#szkielet
import random #do losowania

#podstawowa klasa reprezentujaca pore primerow
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
        

#klasa w której jest algorytm genetyczny
class PrimerDesignGA:
    def __init__(self, dna_sequence, beg_true, end_true, population_size):
        self.dna_sequence = dna_sequence #sekwencja
        self.beg_true = beg_true #poczatek sekwencji do klonowania
        print(self.beg_true)
        self.end_true = end_true #koniec sekwencji do klonowania
        print(self.end_true)
        self.population_size = population_size #wielkosc populacji tj. ilosc primerow
        self.population = self.initialize_population()
        self.display_population() #wyswietla populacje - krok check in!
        
   # inicjalizacja populacji - dodac funkcje sprawdzajaca czy primery sie nie powtarzaja
    def initialize_population(self):
        population = []
        while len(population) < self.population_size:
            fs = random.randint(0, self.beg_true)
            alpha = random.randint(min_primer_length, max_primer_length)
            gamma = random.randint(min_primer_length, max_primer_length)
            beta = random.randint(self.end_true - (fs + alpha), len(self.dna_sequence) - gamma)
            primer_pair = PrimerPair(fs, alpha, beta, gamma)
            if not self.primer_pair_exists(population, primer_pair):
                population.append(primer_pair)
        return population  
        
    def primer_pair_exists(self, population, primer_pair):
        return any(p.fs == primer_pair.fs and p.alpha == primer_pair.alpha and 
                   p.beta == primer_pair.beta and p.gamma == primer_pair.gamma
                   for p in population)
                   
    def display_population(self):
        print("Displaying the entire population of primer pairs:")
        print(f"{'Index':>5} | {'Fs':>5} | {'Fe':>5} | {'Rs':>5} | {'Re':>5} | {'Vector (Fs, Alpha, Beta, Gamma)':>30}")
        print("-" * 70)
        for index, primer_pair in enumerate(self.population):
            print(f"{index:5} | {primer_pair.fs:5} | {primer_pair.fe:5} | {primer_pair.rs:5} | {primer_pair.re:5} | "
                  f"({primer_pair.fs}, {primer_pair.alpha}, {primer_pair.beta}, {primer_pair.gamma})")

    
min_primer_length = 18
max_primer_length = 26

#funkcja do szukania poczatkowego i koncowego indeksu sekwencji do klonowania
def find_target_sequence(dna_sequence, target_sequence):
    start_index = dna_sequence.find(target_sequence)
    if start_index != -1:
        end_index = start_index + len(target_sequence) - 1
        return start_index, end_index
    else:
        return None

whole_dna_sequence = "ATCGTGACTGATCGTACGTACGTAGCTAGTCTAGTCTA"
target_sequence_to_replicate = "GTACGTAGC"
position = find_target_sequence(whole_dna_sequence, target_sequence_to_replicate)
ga = PrimerDesignGA(whole_dna_sequence, position[0], position[1], population_size=50)