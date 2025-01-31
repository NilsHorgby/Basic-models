using Plots
using Random
using StatsBase
using Distributions: Normal
using Unzip

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


#An updated version of the migration function, this time it is non-recursive
#now it also applies selection before the migration function finalizes, improving performace
#The older version is still present below, but it will overflow the stack for large populations
function migration_and_selection(population::Vector{Patch})::Vector{Patch}
    function migration_one_sex(population::Vector{Patch}, sex::Symbol)::Vector{Vector{Individual}}
        current_individuals= [[(ind,patch.location) for ind in getfield(patch, sex)] for patch in population]
        current_individuals = reduce(vcat,reduce(hcat, current_individuals))

        function migrate_one_individual(individual::Tuple{Individual,Int64})::Tuple{Individual,Int64}
            signed_distance::Int64 = Int(round(rand(Normal(0,1.5),1)[1]))
            location = individual[2]
            location = (1 <= location + signed_distance <= number_of_patches) ? location + signed_distance : location
            return (individual[1],location)
        end

        next_individuals = migrate_one_individual.(current_individuals) #we might want this as an Array instead

        #if the the selection is applied here on the Vector{Tuple{Individual, Int64}} we can avoid iterating over all the individuals
        #calculate the fitness of the individuals and return (individual,fitness)
        function calculate_fitness(individual::Tuple{Individual,Int64})::Tuple{Individual,Int64,Float64}
            location = individual[2]
            deliterios_alleles = (location < number_of_patches/2) ? sum(individual[1].genome) : sum(broadcast(~,individual[1].genome))
            return (individual[1],individual[2],1 - selection_cofficient*deliterios_alleles)
        end
        #take a sample of the individuals based on their fitness
        function selection(population::Vector{Tuple{Individual,Int64,Float64}})::Vector{Tuple{Individual,Int64}}
            individuals, locations, fitnesses = unzip(population)
            selected_individuals = sample(zip(individuals,locations)|>collect, Weights(fitnesses),Int(floor(individuals_left_after_selection)), replace=false)
            return selected_individuals
        end

        selected_individuals = selection(calculate_fitness.(next_individuals))
        next_patches = [[] for i in 1:number_of_patches]

        
        for (i,individual) in enumerate(selected_individuals)
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
population = create_population() |> global_mating

migration_and_selection(population)



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
