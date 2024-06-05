'''
This code implements algorithm based on https://academic.oup.com/bioinformatics/article/20/11/1710/300143
User can choose the minimal temperature and maximpal temperature in Celsius degrees for the primers in the code: class PrimerDesignGA
Before using the code,one can specify in driver.py what other parameters will be used throughout the optimization
Authors: Olga Wieromiejczyk, Anna Krzywiecka

'''
import random 
import subprocess
import os
import argparse
#basic class representing a pair of primers
class PrimerPair:
    def __init__(self, fs, alpha, beta, gamma):
        self.fs = fs  #beggining of the forward primer as a number in the dna_sequence
        self.alpha = alpha #length of the forward primer
        self.beta = beta #length of the amplifies sequence between primers
        self.gamma = gamma #length of the reverse primer
        self.fe = fs + alpha
        self.rs = self.fe + beta
        self.re = self.rs + gamma
        self.fitness = None
        self.GC = None
        self.Tmd = None
        self.uni = 0
        self.lengd = None
        self.leng = None
        self.PC = None
        self.Term = None
        self.Sc = None
    def __str__(self):
        return (f'Fs = {self.fs}  ,alpha = {self.alpha}, beta ={self.beta}, gamma = {self.gamma}, self.GC = {self.GC}, self.Tmd = {self.Tmd}, self.uni = {self.uni}, self.lengd = {self.lengd}, self.leng = {self.leng}, self.PC = {self.PC}, self.Term = {self.Term}, self.Sc = {self.Sc}')
    def FITNESS_counting(self):
        '''this is a function which is smaller for less fitting primers and bigger for the best ones. It is the reverse of weight function.'''
        self.fitness = 1/(self.leng  + 3*self.lengd + 3*self.Tmd + 3*self.GC + 3*self.Term + 50*self.uni + 10*self.Sc + 10*self.PC)


#Genetic Algorithm class
class PrimerDesignGA:
    def __init__(self, dna_sequence, beg_true, end_true, population_size, mating_pool, Pe, Pm, max_gen):
        self.dna_sequence = dna_sequence #sequence in which the pair of primers should be found
        self.beg_true = beg_true #start of the target sequence
        self.end_true = end_true #end of the target sequence
        self.population_size = population_size # population size: the amount of pairs of primers in each generations
        self.mating_pool = mating_pool # how many pairs of primers are creating a new population
        self.Pe = Pe # crossing-over likelihood
        self.Pm = Pm # mutation likelihood
        self.max_gen = max_gen # maximal amount of generations
        self.maxtemp = 70
        self.mintemp = 50
        #self.gather_input_info()

        # using the same initial population in gathering the data from different parameters so they can be properly compared
        if os.path.exists("initial_population.txt"):
            self.population = self.read_population_from_file("initial_population.txt")
        else:
            self.population = self.initialize_population()
        self.specifity(0)
        self.population.sort(key=lambda pair: pair.fitness,reverse=True)
        self.new_gen = []
        self.GA()

    def gather_input_info(self):
        '''can be turned on so in each optimalization user can specify minimal and maximal melting temperature of primers in the console '''
        self.mintemp = int(input("Please input the minimum melting temperature: "))
        self.maxtemp = int(input("Please input the maximum melting temperature: "))


    def initialize_population(self):
        population = []
        '''creating a new initial population <<only when this file doesnt exist, because if it does exist, next method is initializing population from file>>'''
        with open("initial_population.txt", "w") as file:
            while len(population) < self.population_size:
                fs = random.randint(0, self.beg_true)
                alpha = random.randint(min_primer_length, max_primer_length)
                gamma = random.randint(min_primer_length, max_primer_length)
                if (len(self.dna_sequence) - gamma - (fs + alpha)) < self.end_true - (fs + alpha):
                    print('ERROR: Invalid primer pair generated')
                else:
                    beta = random.randint(self.end_true - (fs + alpha), len(self.dna_sequence) - gamma - (fs + alpha))
                    primer_pair = PrimerPair(fs, alpha, beta, gamma)
                    if not self.primer_pair_exists(population, primer_pair):
                        self.properties(primer_pair)
                        population.append(primer_pair)
                        file.write(f"{fs},{alpha},{beta},{gamma}\n")
        return population

    def read_population_from_file(self, filename):
        '''loading initial population from file '''
        population = []
        with open(filename, "r") as file:
            for line in file:
                fs, alpha, beta, gamma = map(int, line.strip().split(","))
                primer_pair = PrimerPair(fs, alpha, beta, gamma)
                self.properties(primer_pair)
                population.append(primer_pair)
        return population

    def primer_pair_exists(self, population, primer_pair):
        '''function checking whether a primer pair already exists in population'''
        return any(p.fs == primer_pair.fs and p.alpha == primer_pair.alpha and
                   p.beta == primer_pair.beta and p.gamma == primer_pair.gamma
                   for p in population)

    def display_population(self):
        '''
        way  of monitoring how the population looks like,
        Fs, Fe, Rs, Re are the numbers of nucleotides starting and ending of the forward primer
        and starting and ending of the reverse primer
        '''
        print("Displaying the entire population of primer pairs:")
        print(f"{'Index':>5} | {'Fs':>5} | {'Fe':>5} | {'Rs':>5} | {'Re':>5} | {'Vector (Fs, Alpha, Beta, Gamma)':>30} | {'Fitness':>10}")
        print("-" * 70)
        for index, primer_pair in enumerate(self.population):
            print(f"{index:5} | {primer_pair.fs:5} | {primer_pair.fe:5} | {primer_pair.rs:5} | {primer_pair.re:5} | "
                  f"({primer_pair.fs}, {primer_pair.alpha}, {primer_pair.beta}, {primer_pair.gamma}) | {primer_pair.fitness}")

    def crossover(self, parent1, parent2):
        '''function that creates offspring of given two parent Primer Pairs'''
        R = random.randint(0,15)

        binary_mask = f"{R:04b}"  #to binary

        # crossover - randomly mixing the features Fs, alpha, beta and gamma of parents
        new_fs1, new_fs2 = (parent2.fs if binary_mask[0] == '1' else parent1.fs,
                            parent1.fs if binary_mask[0] == '1' else parent2.fs)
        new_alpha1, new_alpha2 = (parent2.alpha if binary_mask[1] == '1' else parent1.alpha,
                                 parent1.alpha if binary_mask[1] == '1' else parent2.alpha)
        new_beta1, new_beta2 = (parent2.beta if binary_mask[2] == '1' else parent1.beta,
                                parent1.beta if binary_mask[2] == '1' else parent2.beta)
        new_gamma1, new_gamma2 = (parent2.gamma if binary_mask[3] == '1' else parent1.gamma,
                                  parent1.gamma if binary_mask[3] == '1' else parent2.gamma)

        # 2 children with mixed features form parents
        offspring1 = PrimerPair(new_fs1, new_alpha1, new_beta1, new_gamma1)
        offspring2 = PrimerPair(new_fs2, new_alpha2, new_beta2, new_gamma2)
        #checking if they fit 
        if(offspring1.fs + offspring1.alpha + offspring1.beta + offspring1.gamma <= len(self.dna_sequence)):
            if not self.primer_pair_exists(self.population,offspring1):
                if not self.primer_pair_exists(self.new_gen, offspring1):
                    self.properties(offspring1)
                    self.new_gen.append(offspring1)

        if(offspring2.fs + offspring2.alpha + offspring2.beta + offspring2.gamma <= len(self.dna_sequence)):

            if not self.primer_pair_exists(self.population,offspring2):
                if not self.primer_pair_exists(self.new_gen, offspring2):
                    self.properties(offspring2)
                    self.new_gen.append(offspring2)



    def mutate(self, individual):
        '''function that creates offspring from one PrimerPair using mutation'''

        #choosing which parameter to mutate (Fs, alpha, beta, gamma)
        component_to_mutate = random.randint(0, 3)

        if component_to_mutate == 0:
            mutation_value = random.randint(0, self.beg_true)
            mutated_individual = PrimerPair(mutation_value, individual.alpha, individual.beta, individual.gamma)
        elif component_to_mutate == 1:
            mutation_value = random.randint(min_primer_length, max_primer_length)
            mutated_individual = PrimerPair(individual.fs, mutation_value, individual.beta, individual.gamma)
        elif component_to_mutate == 2:
            mutation_value = random.randint(self.end_true - (individual.fs + individual.alpha), len(self.dna_sequence) - individual.gamma - (individual.fs + individual.alpha))
            mutated_individual = PrimerPair(individual.fs, individual.alpha, mutation_value, individual.gamma)
        else:
            mutation_value = random.randint(min_primer_length, max_primer_length)
            mutated_individual = PrimerPair(individual.fs, individual.alpha, individual.beta, mutation_value)
        #makin g sure it is a primer pair that will amplify the target sequence and is in the chosen DNA sequence <which is cut out of genome sequence in function extract_sequence
        if(mutated_individual.fs + mutated_individual.alpha + mutated_individual.beta + mutated_individual.gamma <= len(self.dna_sequence)):
            if not self.primer_pair_exists(self.population, mutated_individual):
                if not self.primer_pair_exists(self.new_gen, mutated_individual):
                    self.properties(mutated_individual)
                    self.new_gen.append(mutated_individual) #adding it to the next population

    def combine_and_sort(self):
        '''this function sorts the primers of 'old' and 'new' population to choose a population of primer pairs with the highest fitness scores'''
        self.specifity(1)
        self.population.extend(self.new_gen)
        self.new_gen.clear()
        self.population.sort(key=lambda pair: pair.fitness,reverse=True)
        return self.population[:min(self.population_size, len(self.population))]

    def new_generation(self):
        ''' function that supervises the process of creating new population from an old one '''
        #Pm is the likelihood of a mutation process, Pe is the likelihood of a crossover process <<functions are above>>
        while len(self.new_gen) < self.mating_pool:
            if(random.random() < self.Pe):
                pair1,pair2 = self.roulette()
                if pair1 is not None and pair2 is not None:
                    self.crossover(pair1,pair2)

            if(random.random() < self.Pm):
                rand_pair = self.population[random.randint(0, self.population_size - 1)]
                self.mutate(rand_pair)

        self.population = self.combine_and_sort()

    def roulette(self):
        '''
        Before crossover, two primer pairs need to be taken of the population to become parents of new primer pairs.
        It is done in a way that favorizes parents with higher fitness scores
        cumulative_score makes a spectrum from0 to 1 where bigger chunks are taken by more fitting individuals.
        Next, a random number from the spectrum is generated.
        Whichever primer pair;s cumulative_score covers this number is selected to mating in crossover above. It happens twice because two parents are needed
        The process guaranttes that the primer pairs with higher fiitness scores are more likely to become parents, because they have a bigget cumulative_score
        '''

        pair_no1 = None
        pair_no2 = None
        if len(self.population)>=2:
            sum_of_fitness = 0
            cumulative_score = 0
            rand1 = random.random()
            rand2 = random.random()
            done = 0

            for pair in self.population:
                sum_of_fitness += pair.fitness

            if sum_of_fitness != 0:
                for pair in self.population:
                    cumulative_score += pair.fitness/sum_of_fitness
                    if pair_no1 is None and rand1 <= cumulative_score:
                        if pair != pair_no2:
                            pair_no1 = pair
                            done +=1

                    if pair_no2 is None and rand2 <= cumulative_score:
                        if pair != pair_no1:
                            pair_no2 = pair
                            done +=1

                    if done == 2:
                        break
                return pair_no1, pair_no2
            else:
                return self.population[random.randint(0, self.population_size - 1)], self.population[random.randint(0, self.population_size - 1)]
        else:
            return None, None


    def GA(self):
        '''
        genetic algorithm is working on a number of populations which are specified in driver.py.
        Each iteration is saving the highest fitness score for a population to a file consisting of parameters chosen in driver.py
        '''
        i = 0
        while i < self.max_gen:
            print(self.population[0].fitness)
            file_name = f"fitness_Pm_{self.Pm}_Pe_{self.Pe}.txt"
            with open(file_name, 'a') as file:
                file.write(f"{self.population[0].fitness}\n")
            self.new_generation()
            i = i + 1





    def properties(self, pair):
        '''
        A function needed to obtain the fitness score of a primer pair,
        it consists of subfunctions counting features of primers specified in the publiation with minor changes,
        that were discussed with a person working with PCR on daily basis and checked in newest publications'''

        seqF = str(self.dna_sequence[pair.fs : pair.fs + pair.alpha])
        #Creating the complementary and reversed sequence of Reverse Primer. This way I receive it's sequence 5'-> 3'
        preseqR = str(self.dna_sequence[pair.fs+pair.alpha +pair.beta : pair.fs+pair.alpha +pair.beta + pair.gamma])
        seqR = str(self.complementary(preseqR))

        #counting the GC content which should be between 40 and 60%

        if (len(seqF) > 0 and len(seqR) > 0):
            FGC = (seqF.count('G')+seqF.count('C'))/len(seqF)
            RGC = (seqR.count('G')+seqR.count('C'))/len(seqR)
            if(FGC >= 0.4 and FGC <= 0.6 and RGC >= 0.4 and RGC <= 0.6):
                pair.GC = 0
            else:
                pair.GC = 1
        else:
            pair.GC = 1


        #Counting the melting temperature difference beatween primers in pair
        #on the basis of https://www.rosalind.bio/en/knowledge/what-formula-is-used-to-calculate-tm 
        #it should be less than 5 Celsius degrees and be between the minimal and maximal temperature specified in the code by user (default: 50, 70)
        if len(seqF)<=13:
            FTM = (seqF.count('G')+seqF.count('C'))*4 + (seqF.count('A')+seqF.count('T'))*2
        else:
            FTM = 64.9+41*(seqF.count('G')+seqF.count('C')-16.4)/(len(seqF))
        if len(seqR)<=13:
            RTM = (seqR.count('G')+seqR.count('C'))*4 + (seqR.count('A')+seqR.count('T'))*2
        else:
            RTM = 64.9+41*(seqR.count('G')+seqR.count('C')-16.4)/(len(seqR))


        if(abs(FTM - RTM) <=5):
            pair.Tmd = 0
            if(FTM <= self.mintemp or FTM >= self.maxtemp or RTM <= self.mintemp or RTM >= self.maxtemp ):
                pair.Tmd = 1
        else:
            pair.Tmd = 1

        #checking the termination which should ba a G or a C for both primers. It can consist of two Gs or Cs but not more
        if len(seqF)>=1:
            if(seqF[-1] in ['G','C']):
                if(len(seqF)>=2):
                    if(seqF[-2] in ['G','C']):
                        if(len(seqF)>=3):
                            if(seqF[-3] not in ['G','C']):
                                pair.Term = 0
                            else:
                                pair.Term = 1
                    else:
                        pair.Term = 0
            else:
                pair.Term = 1
        else:
            pair.Term = 1
        if(len(seqR)>=1):
            if(seqR[-1] in ['G','C']):
                if(len(seqR)>=2):
                    if(seqR[-2] in ['G','C']):
                        if(len(seqR)>=3):
                            if(seqR[-3] in ['G','C']):
                                pair.Term += 1
            else:
                pair.Term += 1
        else:
            pair.Term += 1


        #counting the length difference betaween primers in pair, it is a scale, but cannot be more than 5nn different
        if abs(len(seqF) -len(seqR)) == 5:
            pair.lengd = 0.75
        elif abs(len(seqF) -len(seqR)) < 5 and abs(len(seqF) -len(seqR)) >= 3 :
            pair.lengd = 0.5
        elif abs(len(seqF) -len(seqR)) < 3 and abs(len(seqF) -len(seqR)) > 0 :
            pair.lengd = 0.25
        elif abs(len(seqF) -len(seqR)) == 0:
            pair.lengd = 0
        else:
            pair.lengd = 1
        #analyzing the length of primers. They should be between minimal and maximal primer length, default(18,30)
        if len(seqF) >= min_primer_length and len(seqF) <= max_primer_length and len(seqR) >= min_primer_length and len(seqR) <= max_primer_length:
            pair.leng = 0
        else:
            pair.leng = 1


 
        #checking the self-complementarity
        #hairpins
        min_loop_size = 3
        min_stem_size = 4
        pair.Sc = 0

        if(len(seqF) - 2*min_stem_size - min_loop_size >= 0):
            for i in range(len(seqF)-min_loop_size - 2*min_stem_size+1):
                if self.complementarity_check(seqF[:min_stem_size+i], seqF[min_stem_size + i + min_loop_size:]):
                    pair.Sc = 1
                    break
        if pair.Sc == 0:
            if(len(seqR) - 2*min_stem_size - min_loop_size >= 0):
                for i in range(len(seqR)-min_loop_size - 2*min_stem_size+1):
                    if self.complementarity_check(seqR[:min_stem_size+i], seqR[min_stem_size + i + min_loop_size:]):
                        pair.Sc = 1
                        break

        #linear complementarity

        if(self.complementarity_check(seqF, seqF)):
            pair.Sc = 1
        if(self.complementarity_check(seqR, seqR)):
            pair.Sc = 1


        # checking whether two primers hybrydize together. The minimal number of compatible nucleotides can be changed in function complementarity_check(how_manyTA, hohw_manyCG)

        if(self.complementarity_check(seqF, seqR)):
            pair.PC = 1
        else:
            pair.PC = 0



    @staticmethod
    def complementary(sequence):
        '''function returning an inversed, complementary sequence'''
        compl_sequence =''
        dic = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        for i in sequence:
            compl_sequence += dic[i]
        compl_sequence = compl_sequence[::-1]
        return compl_sequence

    def complementarity_check(self,seq1,seq2):
        '''function returining True or False. It checks whether two sequences are complementary and would hybridize'''
        #deciding how many TA pairs and how many CG pairs will result in a secondary structure
        how_many_TA = 6
        how_many_CG = 4
        summary_how_many = (how_many_TA+how_many_CG)/2

        #artificial change for comparing sequences
        seq2 = self.complementary(seq2)

        if len(seq1) > len(seq2):
            longer = seq1
            shorter = seq2
        else:
            longer = seq2
            shorter = seq1

        for i in range(len(longer)):
            counter_TA = 0
            counter_CG = 0
            for j in range(len(shorter)):
                if i+j >= (len(longer)):
                        break
                if longer[i+j] == shorter[j] and longer[i+j] in ('A', 'T'):
                        counter_TA += 1
                elif longer[i+j] == shorter[j] and longer[i+j] in ('G', 'C'):
                        counter_CG += 1
            if(len(longer) - i-1) < min(how_many_TA, how_many_CG):
                break
            #normalization
            counter_CG = counter_CG*how_many_TA/how_many_CG
            counter_TA = counter_TA*how_many_CG/how_many_TA
            if (counter_CG + counter_TA) >= summary_how_many:
                return True
        return False

    def write_primers_to_fasta(self, sequence, primers, output_file):
        """
        Extracts primer sequences based on their coordinates and writes them to a FASTA file.

        """
        with open(output_file, 'w') as fasta:
            for idx, primer_pair in enumerate(primers):
            # Extracting the forward and reverse primers from the sequence
                fwd_primer = sequence[primer_pair.fs : primer_pair.fs + primer_pair.alpha]  # Convert to zero-based indexing
                rev_primer = sequence[primer_pair.fs + primer_pair.alpha + primer_pair.beta : primer_pair.fs + primer_pair.alpha + primer_pair.beta + primer_pair.gamma]  # Convert to zero-based indexing

            # Writing to FASTA format
                fasta.write(f">{idx}_f\n{fwd_primer}\n")
                #print(f">{idx}_f\n{fwd_primer}\n")
                fasta.write(f">{idx}_r\n{rev_primer}\n")
                #print(f">{idx}_r\n{rev_primer}\n")

    def blast_search(self, fasta_file, blast_db, output_file):
        """
        Runs BLASTN against a given database for each primer in the fasta file.

        Args:
        fasta_file (str): The path to the FASTA file containing primer sequences.
        blast_db (str): The path to the BLAST database.
        output_file (str): The path to save the BLAST output.
        """
        try:
            # Set up the BLASTN command
            command = [
                'blastn',
                '-query', fasta_file,
                '-db', blast_db,
                '-out', output_file,
                '-outfmt', '6',  # Tabular format
                '-perc_identity', '90',  # Minimum percentage identity for matches
                '-qcov_hsp_perc', '90',  # Query coverage per HSP
                '-num_threads', '4',
                '-task', 'blastn-short'  # Optimized for short sequences like primers

            ]

            # Execute the command
            subprocess.run(command, check=True)
            #print("BLASTN search completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"BLASTN search failed: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")

    def count_alignments(self, blast_output_file):
        counts = {}
        with open(blast_output_file, 'r') as file:
            for line in file:
                if line.strip():
                    parts = line.split()
                    primer_name = parts[0]
                    if primer_name in counts:
                        counts[primer_name] += 1
                    else:
                        counts[primer_name] = 1

        return counts

    def write_counts_to_file(self, counts, output_file):
        with open(output_file, 'w') as file:
            for primer, count in counts.items():
                file.write(f"{primer}\t{count}\n")

    def specifity(self, which):
        if which == 0:
            self.write_primers_to_fasta(self.dna_sequence, self.population, "primers_fasta.fasta")
        else:
            self.write_primers_to_fasta(self.dna_sequence, self.new_gen, "primers_fasta.fasta")
        self.blast_search("primers_fasta.fasta","human_genome_db","primers_results.txt")
        results = self.count_alignments("primers_results.txt")
        self.write_counts_to_file(results, "primers_counts.txt")
        os.remove("primers_fasta.fasta")
        os.remove("primers_results.txt")
        os.remove("primers_counts.txt")
        results = list(results.values())
        if which == 0:
            for index,primer in enumerate(self.population):
                if results[index] > 1:
                    primer.uni += 1
                elif results[index] == 0:
                    primer.uni += 1
                if results[index + 1] > 1:
                    primer.uni += 1
                elif results[index + 1] == 0:
                    primer.uni += 1
                primer.FITNESS_counting()
                #print(primer.uni)
        else:
            for index,primer in enumerate(self.new_gen):
                if results[index] > 1:
                    primer.uni += 1
                elif results[index] == 0:
                    primer.uni += 1
                if results[index + 1] > 1:
                    primer.uni += 1
                elif results[index + 1] == 0:
                    primer.uni += 1
                primer.FITNESS_counting()
                #print(primer.uni)

def run_blat(query):
    # Run BLAT
    blat_command = ["blat", "hg38.2bit", query, "blat_output.psl"]
    subprocess.run(blat_command, check=True)

    # Parse the BLAT output to get the best hit coordinates
    best_hit = None
    with open("blat_output.psl", "r") as file:
        lines = file.readlines()
        if len(lines) > 5:  # Skipping header lines in PSL file
            best_hit = lines[5].strip().split()

    # Clean up temporary files
    os.remove("blat_output.psl")
    #for index,field in enumerate(best_hit):
        #print(f"Index {index}: {field}")
    t_starts = list(map(int, best_hit[20].strip(',').split(',')))
    q_sizes = list(map(int, best_hit[18].strip(',').split(',')))
    return best_hit[13], t_starts[0], t_starts[-1] + q_sizes[-1]  # Start and end coordinates of the best hit

def fasta_to_string(file_name):
    try:
        sequence = ''
        with open(file_name, 'r') as file:  # Notice correction in the variable name here
            for line in file:
                line = line.strip()  # Remove trailing newline characters
                if line.startswith('>'):
                    if sequence:  # This checks if sequence is non-empty
                        break  # Stops reading if another header is found (for multi-sequence FASTA files)
                    else:
                        continue  # Skip the header line
                sequence += line  # Append line to sequence
        return sequence  # This line should be indented to be part of the with block
    except Exception as e:  # It's good practice to handle exceptions to know what went wrong
        print("Failed to read file due to:", e)
        return None

def extract_sequence(genome, chrom, start, end):
    # Adjust coordinates to include 500 nucleotides on each side
    start = max(0, start - 500)
    end = end + 500
    # Load the genome using twoBitToFa (if not already in a parseable format)
    subprocess.run(["twoBitToFa", f"-seq={chrom}", f"-start={start}", f"-end={end}", "hg38.2bit", "hg38_cut.fa"], check=True)
    # Extract the sequence
    extended_seq = fasta_to_string("hg38_cut.fa")
    return extended_seq


min_primer_length = 18
max_primer_length = 30
 
def find_target_sequence(dna_sequence, target_sequence):
    '''function to find the first and last index of target sequence'''
    start_index = dna_sequence.find(target_sequence)
    if start_index != -1:
        end_index = start_index + len(target_sequence) - 1
        return start_index, end_index
    else:
        print('The target sequence doesnt exist in the given DNA sequence')
        return None


def parse_args():
    parser = argparse.ArgumentParser(description='Run Primer Design GA with specified Pe and Pm values.')
    parser.add_argument('--Pe', type=float, required=True, help='Crossover probability (Pe)')
    parser.add_argument('--Pm', type=float, required=True, help='Mutation probability (Pm)')
    return parser.parse_args()


def main():
    '''Main function to run the GA with specified parameters'''
    args = parse_args()

    if os.path.exists("hg38_cut.fa"):
        extended_sequence = fasta_to_string("hg38_cut.fa")
    else:
        target_sequence = "target.fasta"
        chrom, start, end = run_blat(target_sequence)
        extended_sequence = extract_sequence("hg38.fa", chrom, start, end)

    extended_sequence = extended_sequence.upper()

    # Initialize and run the GA
    ga = PrimerDesignGA(
        extended_sequence,
        1000,
        len(extended_sequence) - 1000,
        population_size=200,
        mating_pool=80,
        Pe=args.Pe,
        Pm=args.Pm,
        max_gen=100
    )


if __name__ == "__main__":
    main()