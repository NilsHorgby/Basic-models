# Explenation 
# This is a diploid model with sex. 

# The allels are not currently unique, 
# an "A" allel can end up as an "a" allel in the next generation
# "B" or "b" will never end up as "A" or "a" and wise versa

# The state of an individuals genetic make up is designated with an int from 0 - 3
# This is a binary representation of both allels. 
# 0 : Bb, 1 : Ba, 2 : Ab, 3 : Aa



#population of males : Vector{Int}
#   size N, values from 0-3
#population of females : Vector{Int}
#   size N, values from 0-3

#take a random selection of males and females and pair them up randomly 
#take a random selaction of one allel from each first_parent
#make the child female or male randomly

using StatsBase
N_male::Int = 100
N_female::Int = 100
p::Float64 = 0.55 #intitial Frequency of Allel a 
total_time::Int = 1000
sampling_frequency::Int = 1


#    |A      |B      
#A   |p²     |p*(1-p)
#B   |p*(1-p)|(1-p)²

function make_first_gen(p::Float64, population_size::Int)
    weights = Weights([(1-p)^2,(1-p)*p,(1-p)*p,p^2])
    population::Vector{Int8} = sample([0,1,2,3],weights,population_size,replace = true )
    return population
end

#function takes one parent, this function is broadcasted over the whole population
#return one of the parent's allels at random
function make_sex_cells(parent)
    if parent == 0
        return 0
    elseif parent <= 2
        return rand([0,1])
    else
        return 1
    end
end


function make_new_generation(female_population::Vector{Int8},male_population::Vector{Int8}, population_size::Int)
    #Designate parents. One individual can be a parent multiple times
    female_parents::Vector{Int8} = sample(female_population, 2*population_size)
    male_parents::Vector{Int8} = sample(male_population, 2*population_size)

    #Create haploid "sex cell" population
    eggs::Vector{Int8} = make_sex_cells.(female_parents)
    sperm::Vector{Int8} = make_sex_cells.(male_parents)

    #Combining the sex cells into the next generation
    #For simplisities sake the females will always decide the large letter Allel
    next_generation::Vector{Int8} = eggs*2 + sperm 

    #Deciding what individual will end upp female/male
    #The total population size will be the same as the previous generation
    #The number of females and males vary randomly. The number of each sex has no effect on the next generation.
    female_mask::BitVector = rand(Bool, length(next_generation))
    male_mask::BitVector = .!female_mask
    
    next_females::Vector{Int8} = next_generation[female_mask]
    next_males::Vector{Int8} = next_generation[male_mask]
    return next_females, next_males
end

function allel_frequency(population)
    p::Float64 = (sum(population .== 1) + sum(population .== 2) + 2*sum(population .== 3))/(2*length(population))
    return p
end

function fisher_wright_single_run(p::Float64, N_female::Int, N_male::Int, total_time::Int, sampling_frequency::Int)
    #intitialize both males and female populations 
    population_size::Int = N_female + N_male
    female_population::Vector{Int8} = make_first_gen(p, N_female)
    male_population::Vector{Int8} = make_first_gen(p, N_male)
    allel_frequency_over_time::Vector{Tuple{Float64,Float64}} = [(p,p)]

    for generation in 1:total_time
        female_population, male_population = make_new_generation(female_population,male_population, population_size)
        if generation % sampling_frequency == 0
            push!(allel_frequency_over_time,allel_frequency.((female_population,male_population)))
        end
    end

    return allel_frequency_over_time
end

fisher_wright_single_run(p,N_female,N_male,total_time,sampling_frequency)