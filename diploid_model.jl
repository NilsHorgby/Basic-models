#Diploid Model

#original model has individuls where each individual is a column with two entries (two homologous chromosomes) 
#possibly it would be better to represent the diffrent possible states in binary/ ints between 0 and 3, 0 : [0,0], 1 : [0,1] ... 3: [1,1]
using StatsBase

population_size = 1000
p = 0.55 #intitial Frequency of Allel a 
total_time = 1
sampling_frequency = 1

function random_parent(this_generation, population_size)
    indecies = rand(1:population_size, population_size)
    return this_generation[:,indecies]
end

function time_step(this_generation, population_size)
    first_parent::Matrix{Bool} = random_parent(this_generation, population_size)
    second_parent = random_parent(this_generation, population_size)
    first_allel = sample.(eachslice(first_parent, dims = 2))
    second_allel = sample.(eachslice(second_parent, dims = 2))
    next_generation = reshape(vcat(first_allel, second_allel), 2, population_size)
    return next_generation
end

function allel_frequency(population,population_size )
    return sum(eachslice(population, dims = 1)[1])/population_size
end

function fisher_wright_single_run(p::Float64, population_size::Int, total_time::Int, sampling_frequency::Int)
    first_allel::Vector{Bool} = sample([0,1],Weights([1-p,p]), population_size, replace = true) #whole populations' allel a-s
    second_allel::Vector{Bool} = zeros(population_size) #whole populations' allel b-s
    population::Matrix{Bool} = reshape(vcat(first_allel, second_allel), 2, population_size) # N by 2 matrix
    allel_frequency_over_time::Vector{Float64} = [p] #intitialize with initial allel frequency

    for generation in 1:total_time
        population = time_step(population, population_size)
        if generation % sampling_frequency == 0
            push!(allel_frequency_over_time, allel_frequency(population, population_size)) #something very strange going on with the frequencies
        end
    end
    
    return population, total_time, allel_frequency_over_time
end
fisher_wright_single_run(p, population_size, total_time, sampling_frequency)


