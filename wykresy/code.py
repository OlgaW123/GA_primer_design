#szkielet
import random #do losowania
import subprocess
import os
import argparse
#podstawowa klasa reprezentujaca pare primerow
class PrimerPair:
    def __init__(self, fs, alpha, beta, gamma):
        self.fs = fs
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.fe = fs + alpha
        self.rs = self.fe + beta
        self.re = self.rs + gamma
        self.fitness = None
        self.GC = None
        self.Tmd = None
        self.uni = 0
        self.lengd = None
        self.leng = None
        #self.R = None
        self.PC = None
        self.Term = None
        self.Sc = None
    def __str__(self):
        return (f'Fs = {self.fs}  ,alpha = {self.alpha}, beta ={self.beta}, gamma = {self.gamma}, self.GC = {self.GC}, self.Tmd = {self.Tmd}, self.uni = {self.uni}, self.lengd = {self.lengd}, self.leng = {self.leng}, self.PC = {self.PC}, self.Term = {self.Term}, self.Sc = {self.Sc}')
    def FITNESS_counting(self):
        self.fitness = 1/(self.leng  + 3*self.lengd + 3*self.Tmd + 3*self.GC + 3*self.Term + 50*self.uni + 10*self.Sc + 10*self.PC)


#klasa w której jest algorytm genetyczny
class PrimerDesignGA:


    def __init__(self, dna_sequence, beg_true, end_true, population_size, mating_pool, Pe, Pm, max_gen):
        self.dna_sequence = dna_sequence #sekwencja
        self.beg_true = beg_true #poczatek sekwencji do klonowania
        #print(self.beg_true)
        self.end_true = end_true #koniec sekwencji do klonowania
        #print(self.end_true)
        self.population_size = population_size #wielkosc populacji tj. ilosc primerow
        self.mating_pool = mating_pool #ilosc potomstwa
        self.Pe = Pe #prawdopodobienstwo crossing-over
        self.Pm = Pm #prawdopodobienstwo mutacji
        self.max_gen = max_gen #maksymalna ilosc generacji
        #self.restriction_sequences = []
        self.maxtemp = 50
        self.mintemp = 70
        #self.gather_input_info()
        if os.path.exists("initial_population.txt"):
            self.population = self.read_population_from_file("initial_population.txt")
        else:
            self.population = self.initialize_population()
        self.specifity(0)
        self.population.sort(key=lambda pair: pair.fitness,reverse=True)
        self.new_gen = []
        self.GA()

    def gather_input_info(self):
        self.mintemp = int(input("Please input the minimum melting temperature: "))
        self.maxtemp = int(input("Please input the maximum melting temperature: "))



   # inicjalizacja populacji - dodac funkcje sprawdzajaca czy primery sie nie powtarzaja   ------ wydaje mi sie ze sprawdzasz przed dodaniem
    def initialize_population(self):
        population = []
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
        population = []
        with open(filename, "r") as file:
            for line in file:
                fs, alpha, beta, gamma = map(int, line.strip().split(","))
                primer_pair = PrimerPair(fs, alpha, beta, gamma)
                self.properties(primer_pair)
                population.append(primer_pair)
        return population

    def primer_pair_exists(self, population, primer_pair):
        return any(p.fs == primer_pair.fs and p.alpha == primer_pair.alpha and
                   p.beta == primer_pair.beta and p.gamma == primer_pair.gamma
                   for p in population)

    def display_population(self):
        print("Displaying the entire population of primer pairs:")
        print(f"{'Index':>5} | {'Fs':>5} | {'Fe':>5} | {'Rs':>5} | {'Re':>5} | {'Vector (Fs, Alpha, Beta, Gamma)':>30} | {'Fitness':>10}")
        print("-" * 70)
        for index, primer_pair in enumerate(self.population):
            print(f"{index:5} | {primer_pair.fs:5} | {primer_pair.fe:5} | {primer_pair.rs:5} | {primer_pair.re:5} | "
                  f"({primer_pair.fs}, {primer_pair.alpha}, {primer_pair.beta}, {primer_pair.gamma}) | {primer_pair.fitness}")

    def crossover(self, parent1, parent2):

        R = random.randint(0,15)

        binary_mask = f"{R:04b}"  #na binarne

        # crossover
        new_fs1, new_fs2 = (parent2.fs if binary_mask[0] == '1' else parent1.fs,
                            parent1.fs if binary_mask[0] == '1' else parent2.fs)
        new_alpha1, new_alpha2 = (parent2.alpha if binary_mask[1] == '1' else parent1.alpha,
                                 parent1.alpha if binary_mask[1] == '1' else parent2.alpha)
        new_beta1, new_beta2 = (parent2.beta if binary_mask[2] == '1' else parent1.beta,
                                parent1.beta if binary_mask[2] == '1' else parent2.beta)
        new_gamma1, new_gamma2 = (parent2.gamma if binary_mask[3] == '1' else parent1.gamma,
                                  parent1.gamma if binary_mask[3] == '1' else parent2.gamma)

        # dzieci
        offspring1 = PrimerPair(new_fs1, new_alpha1, new_beta1, new_gamma1)
        offspring2 = PrimerPair(new_fs2, new_alpha2, new_beta2, new_gamma2)

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

        if(mutated_individual.fs + mutated_individual.alpha + mutated_individual.beta + mutated_individual.gamma <= len(self.dna_sequence)):
            if not self.primer_pair_exists(self.population, mutated_individual):
                if not self.primer_pair_exists(self.new_gen, mutated_individual):
                    self.properties(mutated_individual)
                    self.new_gen.append(mutated_individual)

    def combine_and_sort(self):
        self.specifity(1)
        self.population.extend(self.new_gen)
        self.new_gen.clear()
        self.population.sort(key=lambda pair: pair.fitness,reverse=True)
        return self.population[:min(self.population_size, len(self.population))]

    def new_generation(self):
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
            #print(sum_of_fitness)


            if sum_of_fitness != 0:
                for pair in self.population:
                    #print(pair)
                    cumulative_score += pair.fitness/sum_of_fitness
                    #print(f'pairfitness: {pair.fitness}')
                    #print(f'cumul {cumulative_score}')
                    if pair_no1 is None and rand1 <= cumulative_score:
                        if pair != pair_no2:
                            pair_no1 = pair
                            done +=1
                            #print(f'pair1{pair_no1}')

                    if pair_no2 is None and rand2 <= cumulative_score:
                        if pair != pair_no1:
                            pair_no2 = pair
                            #print(f'pair2{pair_no2}')
                            done +=1
                    if done == 2:
                        break
                return pair_no1, pair_no2
            else:
                #tu w sumie moznaby przerwa algorytm np dodac jakis argument do GA jakby czy np nie jest idealnie i mozna dac tam break - watpie ze to sie wydarzy ale no
                return self.population[random.randint(0, self.population_size - 1)], self.population[random.randint(0, self.population_size - 1)]
        else:
            return None, None


    def GA(self):
        i = 0
        while i < self.max_gen:
            print(self.population[0].fitness)
            file_name = f"fitness_Pm_{self.Pm}_Pe_{self.Pe}.txt"
            with open(file_name, 'a') as file:
                file.write(f"{self.population[0].fitness}\n")
            self.new_generation()
            i = i + 1
            #self.display_population()





    def properties(self, pair):

        seqF = str(self.dna_sequence[pair.fs : pair.fs + pair.alpha])
        #here I create the complementary and reversed sequence of Reverse Primer. This way I receive it's sequence 5'-> 3'
        preseqR = str(self.dna_sequence[pair.fs+pair.alpha +pair.beta : pair.fs+pair.alpha +pair.beta + pair.gamma])
        seqR = str(self.complementary(preseqR))

        #gc_counting
        if (len(seqF) > 0 and len(seqR) > 0):
            FGC = (seqF.count('G')+seqF.count('C'))/len(seqF)
            RGC = (seqR.count('G')+seqR.count('C'))/len(seqR)
            if(FGC >= 0.4 and FGC <= 0.6 and RGC >= 0.4 and RGC <= 0.6):
                pair.GC = 0
            else:
                pair.GC = 1
        else:
            pair.GC = 1


        #Tmd_counting
        #on the basis of https://www.rosalind.bio/en/knowledge/what-formula-is-used-to-calculate-tm
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

        #Term_counting
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


        #lengd_counting
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
        #leng_counting
        if len(seqF) >= min_primer_length and len(seqF) <= max_primer_length and len(seqR) >= min_primer_length and len(seqR) <= max_primer_length:
            pair.leng = 0
        else:
            pair.leng = 1


        #Sc_counting
        #here I check the self-complementarity
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

        #linear

        if(self.complementarity_check(seqF, seqF)):
            pair.Sc = 1
        if(self.complementarity_check(seqR, seqR)):
            pair.Sc = 1


        #PC_counting
        #here we check whether two primers hybrydize together. The minimal number of compatible nucleotides is in the variable min_compl_size

        if(self.complementarity_check(seqF, seqR)):
            pair.PC = 1
        else:
            pair.PC = 0



    @staticmethod
    def complementary(sequence):
        compl_sequence =''
        dic = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        for i in sequence:
            compl_sequence += dic[i]
        compl_sequence = compl_sequence[::-1]
        return compl_sequence

    def complementarity_check(self,seq1,seq2):
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
    # Adjust coordinates to include 1000 nucleotides on each side
    start = max(0, start - 1000)
    #print(start)
    end = end + 1000
    #print(end)
    # Load the genome using twoBitToFa (if not already in a parseable format)
    subprocess.run(["twoBitToFa", f"-seq={chrom}", f"-start={start}", f"-end={end}", "hg38.2bit", "hg38_cut.fa"], check=True)
    # Extract the sequence
    extended_seq = fasta_to_string("hg38_cut.fa")
    return extended_seq


#czemu nie damy tego wczesniej? bo jest to uzywane w funkcjach powyzej
            #MUSIMY ZMIENIC!!!!!!
min_primer_length = 18
max_primer_length = 30
#funkcja do szukania poczatkowego i koncowego indeksu sekwencji do klonowania
def find_target_sequence(dna_sequence, target_sequence):
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

# Main function to run the GA with specified parameters
def main():
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
        population_size=500,
        mating_pool=100,
        Pe=args.Pe,
        Pm=args.Pm,
        max_gen=100
    )


if __name__ == "__main__":
    main()
