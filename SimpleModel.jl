using Plots
using Random
using StatsBase

#Constants
const individuals_per_patch = 100
const number_of_patches = 50
const individuals_born_per_time_step = 1000
const individuals_left_after_selection = 100
const selection_cofficient = 0.1

#An individual has a genome, a sex and an id
struct Individual
    genome:: Tuple{Bool, Bool}
    sex:: Bool #false = female, true = male
    id:: Int
end

#A patch has a location and a list of individuals
struct Patch
    females:: Vector{Individual}
    males:: Vector{Individual}
    location:: Int
end


#create the initial population
function create_individual(id,sex)
    genome = (rand(Bool), rand(Bool))
    return Individual(genome,sex,id)
end

#create the population of a single patch
function create_subpopulation(location)
    females = [create_individual(location*individuals_per_patch/2*1 + i,true) for i in 1:Int(floor(individuals_per_patch/2))]
    males = [create_individual(location*individuals_per_patch/2*2 + i,false) for i in 1:Int(floor(individuals_per_patch/2))]
    return Patch(females,males,location)
end

function create_population()
    patches = [create_subpopulation(i) for i in 1:number_of_patches]
    return patches
end

getfield(create_subpopulation(1),:males)


#this mating function pools all gametes in one pool
function local_mating(patch::Patch)::Patch
    #select a random allele from each parent
    female_gametes = sample((individual -> individual.genome[rand([1,2])]).(patch.females), Int(floor(individuals_born_per_time_step)))
    male_gametes = sample((individual -> individual.genome[rand([1,2])]).(patch.males), Int(floor(individuals_born_per_time_step)))
    
    #create a new individuals with the selected alleles
    offspring_genomes = zip(female_gametes, shuffle(male_gametes))|>collect
    female_genomes = offspring_genomes[1:Int(floor(individuals_born_per_time_step/2))]
    male_genomes = offspring_genomes[Int(floor(individuals_born_per_time_step/2)+1):end]

    #create the offspring, id is unique for each individual in the population
    female_offspring = [Individual(genome, false, patch.location*individuals_born_per_time_step/2*1 + i) for (i, genome) in enumerate(female_genomes)]
    male_offspring = [Individual(genome, true, patch.location*individuals_born_per_time_step/2*2 + i) for (i, genome) in enumerate(male_genomes)]

    #return the new patch
    return Patch(female_offspring,male_offspring, patch.location)
end

#Apply the local mating function to all patches in the population
function global_mating(population::Vector{Patch})::Vector{Patch}
    return local_mating.(population)
end

using Distributions: Normal


struct SexSpecificPatch
    individuals::Vector{Individual}
    sex::Bool
    location::Int
end

#An updated version of the migration function, this time it is non-recursive
#The older version is still present below, but it will overflow the stack for large populations
function migration(population::Vector{Patch})::Vector{Patch}
    function migration_one_sex(population::Vector{Patch}, sex::Symbol)::Vector{Vector{Individual}}
        current_individuals= [[(ind,patch.location) for ind in getfield(patch, sex)] for patch in population]
        current_individuals = reduce(vcat,reduce(hcat, current_individuals))

        function migrate_one_individual(individual::Tuple{Individual,Int64})::Tuple{Individual,Int64}
            signed_distance::Int64 = Int(round(rand(Normal(0,1.5),1)[1]))
            location = individual[2]
            location = (1 <= location + signed_distance <= number_of_patches) ? location + signed_distance : location
            return (individual[1],location)
        end

        next_individuals = migrate_one_individual.(current_individuals)
        next_patches = [[] for i in 1:number_of_patches]

        #if the the selection is applied here on the Vector{Tuple{Individual, Int64}} we can avoid iterating over all the individuals
        for (i,individual) in enumerate(next_individuals)
            push!(next_patches[individual[2]],individual[1])
        end
        
        return next_patches
    end

    female_patches = migration_one_sex(population,:females)
    male_patches = migration_one_sex(population,:males)

    function combine_males_and_females(males::Vector{Individual},females::Vector{Individual},location::Int)::Patch
        return Patch(males,females,location)
    end

    return combine_males_and_females.(female_patches,male_patches,1:50)
end



#this function as the side effect of making the current population empty, the population must be assigned to the next population
#This function is a bit of recersive mess, but it should work
function recursive_migration(population::Vector{Patch})::Vector{Patch}
    current_females = [SexSpecificPatch(patch.females,false, patch.location) for patch in population]
    current_males = [SexSpecificPatch(patch.males,true, patch.location) for patch in population]
    next_females::Vector{SexSpecificPatch} = [SexSpecificPatch([],false,population[i].location)
                                              for i=eachindex(population)]
    next_males::Vector{SexSpecificPatch} = [SexSpecificPatch([],true,population[i].location)
                                            for i=eachindex(population)]
    #migration one by one:

    function put_within_boundaries(location,signed_distance)
        new_patch::Int64 = location + signed_distance
        if new_patch <= 0
            return 1
        elseif new_patch > number_of_patches
            return number_of_patches
        else 
            return new_patch
        end
    end
    
    function migrate_one_individual(population::Vector{SexSpecificPatch},next_population::Vector{SexSpecificPatch})
        individual::Individual= population[1].individuals[1]
        signed_distance::Int64 = Int(round(rand(Normal(0,1.5),1)[1]))
        deleteat!(population[1].individuals,1)
        
        #Boundary condition: if the individual tries to go outside the area, it doesn't move
        new_patch = put_within_boundaries(population[1].location,signed_distance)
        append!(next_population[new_patch].individuals,[individual])  
        if length(population[1].individuals) == 0
            deleteat!(population,1)
        end
        return population,next_population
    end
    
    #we should be able to do tail recursion here
    function migrate(population::Vector{SexSpecificPatch},next_population::Vector{SexSpecificPatch})
        if population == []
            return next_population
        else
            return migrate(migrate_one_individual(population,next_population)...)
        end
    end

    
    function combine_males_and_females(females::SexSpecificPatch,males::SexSpecificPatch)::Patch
        @assert females.location == males.location "The locations are not the same"
        return Patch(females.individuals,males.individuals, females.location)
    end
    return combine_males_and_females.(migrate(current_females,next_females),migrate(current_males,next_males))
end

#Selection

function calculate_fitness(individual::Individual, location::Int64)::Float64
    #returns:
        #1 if the individual has no deliterious alleles
        #1 - selection_cofficient if the individual has one deliterious allele
        #1 - 2*selection_cofficient if the individual has two deliterious alleles
    deliterios_alleles = (location < number_of_patches/2) ? sum(individual.genome) : sum(broadcast(~,individual.genome)) #if the individual is in the first half of the patches, the first allele is the deliterious one, otherwise the second allele is the deliterious one
    return 1 - selection_cofficient*deliterios_alleles
end



function selection(population::Vector{Patch})::Vector{Patch}
    function selection_by_patch(patch::Patch)::Patch
        calculate_fitness_partial = ((ind) -> calculate_fitness(ind,patch.location))
        male_fitnesses = calculate_fitness_partial.(patch.males)
        female_fitnesses = calculate_fitness_partial.(patch.females)
        males_remaining = sample(patch.males, Weights(male_fitnesses),Int(floor(individuals_left_after_selection/2)), replace=false)
        female_remaining = sample(patch.females, Weights(female_fitnesses),Int(floor(individuals_left_after_selection/2)), replace=false)
        return Patch(female_remaining,males_remaining,patch.location)
    end
    return selection_by_patch.(population)
end


function get_allele_frequencies(population::Vector{Patch}, generations::Int64)::Vector{Float64}
    for n in 1:100
        population = global_mating(population) |> migration |> selection
    end
    return mean.(((patch) -> ((male) -> sum(male.genome)).(patch.males)).(population))
end

function predicted_allele_frequencie(patch)::Float64

    x = patch - number_of_patches/2
    if x >= 0
        return -1 + 3 * tanh(sqrt(selection_cofficient)*x/(2*1.5) + atanh(sqrt(2/3)))^2
    else
        return 3 - 3 * tanh(-sqrt(selection_cofficient)*x/(2*1.5) + atanh(sqrt(2/3)))^2
    end
end

population = create_population()



frequencies = [get_allele_frequencies(create_population(),100) for i in 1:100]

plot(frequencies,linestyle = :dashdot, color = "gray", label = "")
plot!(mean.(collect(eachrow(reduce(hcat, frequencies)))), 
    title = "Allele frequency in Simple Model, mean of 100 runs",
    xlabel = "Patch",
    ylabel = "Allele frequency",
    line = (3,:black,:dashdot),
    label = "Mean Observed Allele Frequencies")

predicted_frequencies = [predicted_allele_frequencie(i) for i in range(1,Int(number_of_patches),length = 1000)]
plot!(range(1,Int(number_of_patches),length = 1000),predicted_frequencies,
    label = "Predicted allele frequencies",
    line=(3,:red))
savefig("BasicModel_with_pred.png")




#Main function
function main()
    population = create_population()
    for n in 1:100
        population = global_mating(population) |> migration |> selection
    end
    plot(mean.(((patch) -> ((male) -> sum(male.genome)).(patch.males)).(population)),
        title = "Allele frequency in Simple Model",
        xlabel = "Patch",
        ylabel = "Allele frequency")
    savefig("BasicModel.png") 
end
main()
########################################################
#Tests


#Test how far the individuals have migrated in one time step
function test_migration_distances()::Tuple{Vector{Float64},Vector{Float64}}
    population = create_population()
    population = migration(population)
    function exctract_origin(patch::Patch)::Vector{Int64} #only for males 
        ids::Vector{Int64} = ((ind) -> ind.id).(patch.males)
        location_of_origin::Vector{Int64} = ((id)-> (id - mod(id-1,individuals_per_patch/2)-1)/(individuals_per_patch/2)).(ids)
        return location_of_origin
    end

    function migration_distance(population::Vector{Patch})
        return [((origin) -> origin - i).(origins) for (i,origins) in enumerate(exctract_origin.(population))]
    end


    mean_signed_distance::Vector{Float64} = mean.(migration_distance(population))
    
    mean_absolute_distance::Vector{Float64} = mean.(((x)->broadcast(^,x,2)).(migration_distance(population))) #should on average be the variance of the normal distribution

    return mean_signed_distance,mean_absolute_distance
end


mean_signed_distance,mean_absolute_distance = test_migration_distances()
sqrt(mean(mean_absolute_distance[3:45])) #should be ~1.5, the standard deviation of the normal distribution, the first and last patch will have diffrent values, so we ignore them
mean(mean_signed_distance) #should be ~0, the mean of the normal distribution

#Tests for selection

function test_fittness_calculations()
    selection_cofficient = 0.1
    return (calculate_fitness(Individual((false,false),true,1),1) == 1.0) & (calculate_fitness(Individual((false,true),true,1),1) == (1.0 - selection_cofficient)) & (calculate_fitness(Individual((true,true),true,1),1) == (1.0 - 2*selection_cofficient))
end
@assert test_fittness_calculations() "The fitness calculations are not correct"

#Test that the population size is constant

function test_constant_population_size(population::Vector{Patch})::Bool
    function number_of_individuals_in_patch(patch::Patch)::Int
        return length(patch.males) + length(patch.females)
    end

    return sum(number_of_individuals_in_patch.(population)) == individuals_per_patch*number_of_patches 
end
#Test that all patches are present
function test_all_patches_present(population::Vector{Patch})::Bool
    return (patch -> patch.location).(population) == 1:number_of_patches
    
end


population = create_population()
@assert test_constant_population_size(population) == true "The population size is not constant"
population = migration(population)
@assert test_constant_population_size(population) == true "The population size is not constant"
population = global_mating(population)
@assert test_all_patches_present(population) "Not all patches are present"
