import random
from functools import partial
from itertools import cycle, chain, repeat
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import requests
import re
import time

from deap import base
from deap import creator
from deap import tools

mutated =  "ATGGCCAAGAAGACCTCATCTAAGGGCAAGCTGCCCCCCGGCCCCTTCCCCCTGCCCATCATCGGCAACCTGTTCGACTTCGAGCTGAAGAACATCCCCAAGTCCTTCACCCGGCTGGCCCAGCGGTTCTCCCCCGTCTTCACCCTGTATGTCAACTCCCAGCGGATGGTCGTCATGCACGGCTATAAGGCCGTCAAGGAGGCCCTGCTGGACTATAAGGACGAGTTCTCCGGCCGGGGCGACCTGCCCTCCGACCACGCCCACCGGGACCGGGGCATCATCTTCAACAACGGCCCCACCTGGAAGGACATCCGGCGGTTCTCCCTGACCACCCTGCGGAACTATGGCATGGGCAAGCAGGGCAACGAGTCCCGGATCCAGCGGGAGAAGCACTTCCTGCTGGAGGCCCTGCGGAAGACCCAGGGCCAGCCCTTCGACCCCACCTTCCTGATCGGCTGCGCCCCCTGCAACGTCATCGCCGACATCCTGTTCCGGAAGCACTTCGACTATAACGACGAGAAGTTCCTGCGGCTGATGTATCTGTTCAACGAGAACTTCCACCTGCGGTCCACCCCCTGGCTGCAGCTGTATAACAACTTCCCCTCCTTCCTGCACTATCTGCCCGGCTCCCACCGGAAGGTCATCAAGAACGTCGCCGAGGTCAAGGAGTATGTCTCCCAGCGGGTCGAGGAGCACCACCAGTCCCTGGACCCCAACTGCCCCCGGGACCTGACCGACTGCCTGCTGGTCGAGATGGAGAAGGAGAAGCACTCCGCCGAGCGGCTGTATACCATGGACGGCATCACCGTCACCGTCGCCAAGCTGTTCTTCGCCGGCACCGAGACCACCTCCACCACCCTGCGGTATGGCCTGCTGATCCTGATGAAGTATCCCGAGATCGAGGAGAAGGACCACGAGGAGATCGACCGGGTCATCGGCCCCTCCCGGATCCCCGCCATCAAGGACCGGCAGGAGATGCCCTATATGGACGCCGTCGTCCACGAGATCCAGCGGTTCATCACCCTGATGCCCTCCAACTGCCCCCACTTCGCCACCCGGGACACCATCTTCCGGGGCTATCTGATCCCCAAGGAGACCGTCGTCGTCCCCACCCTGGACTCCGTCCTGTATGACAACCAGGAGTTCCCCGACCCCGAGAAGTTCGACCCCGAGGAGTTCCTGAACGAGAACGGCAAGTTCAAGTATTCCGACTATTTCAAGCCCTTCTCCACCGGCAAGCGGGTCTGCGCCGGCGAGGGCCTGGCCCGGATGGAGCTGTTCCTGCTGCTGTGCGCCATCCTGCAGCACTTCAACCTGAAGCCCCTGGTCGACCCCAAGGACATCGACCTGTCCCCCCGGCACATCGGCCACGGCTGCATCCCCCCCCGGTATCGGCTGGAGGTCATCCCCCGGTCCTGA"
original = "ATGGCCAAGAAGACCTCATCTAAGGGCAAGCTGCCCCCCGGCCCCTTCCCCCTGCCCATCATCGGCAACCTGTTCCAGCTGGAGCTGAAGAACATCCCCAAGTCCTTCACCCGGCTGGCCCAGCGGTTCGGCCCCGTCTTCACCCTGTATGTCGGCTCCCAGCGGATGGTCGTCATGCACGGCTATAAGGCCGTCAAGGAGGCCCTGCTGGACTATAAGGACGAGTTCTCCGGCCGGGGCGACCTGCCCGCCTTCCACGCCCACCGGGACCGGGGCATCATCTTCAACAACGGCCCCACCTGGAAGGACATCCGGCGGTTCTCCCTGACCACCCTGCGGAACTATGGCATGGGCAAGCAGGGCAACGAGTCCCGGATCCAGCGGGAGGCCCACTTCCTGCTGGAGGCCCTGCGGAAGACCCAGGGCCAGCCCTTCGACCCCACCTTCCTGATCGGCTGCGCCCCCTGCAACGTCATCGCCGACATCCTGTTCCGGAAGCACTTCGACTATAACGACGAGAAGTTCCTGCGGCTGATGTATCTGTTCAACGAGAACTTCCACCTGCTGTCCACCCCCTGGCTGCAGCTGTATAACAACTTCCCCTCCTTCCTGCACTATCTGCCCGGCTCCCACCGGAAGGTCATCAAGAACGTCGCCGAGGTCAAGGAGTATGTCTCCGAGCGGGTCAAGGAGCACCACCAGTCCCTGGACCCCAACTGCCCCCGGGACCTGACCGACTGCCTGCTGGTCGAGATGGAGAAGGAGAAGCACTCCGCCGAGCGGCTGTATACCATGGACGGCATCACCGTCACCGTCGCCGACCTGTTCTTCGCCGGCACCGAGACCACCTCCACCACCCTGCGGTATGGCCTGCTGATCCTGATGAAGTATCCCGAGATCGAGGAGAAGCTGCACGAGGAGATCGACCGGGTCATCGGCCCCTCCCGGATCCCCGCCATCAAGGACCGGCAGGAGATGCCCTATATGGACGCCGTCGTCCACGAGATCCAGCGGTTCATCACCCTGGTCCCCTCCAACCTGCCCCACGAGGCCACCCGGGACACCATCTTCCGGGGCTATCTGATCCCCAAGGGCACCGTCGTCGTCCCCACCCTGGACTCCGTCCTGTATGACAACCAGGAGTTCCCCGACCCCGAGAAGTTCAAGCCCGAGCACTTCCTGAACGAGAACGGCAAGTTCAAGTATTCCGACTATTTCAAGCCCTTCTCCACCGGCAAGCGGGTCTGCGCCGGCGAGGGCCTGGCCCGGATGGAGCTGTTCCTGCTGCTGTGCGCCATCCTGCAGCACTTCAACCTGAAGCCCCTGGTCGACCCCAAGGACATCGACCTGTCCCCCATCCACATCGGCTTCGGCTGCATCCCCCCCCGGTATAAGCTGTGCGTCATCCCCCGGTCCTGA"
#the_best = "MAKKTSSKGKLPPGPFPLPIIGNLLDFELKNIPKSFTRLAQRMGPVFTKYVNSQRMVVMHGYKAVKEALLDYKDEFSGRGDLPSDHAHRDRGIIFNNGPTWKDIRRFSLTTLRNYGMGKQGNESRIQREKHFLLEALRKTQGQPFDPTFLIGCAPCNVIADILFRKHFDYNDEKFLRLMYLFNENFHLRSTPWLQLYNNFPSFLHYLPGSHRKVIKNVAEVKEYVSQRVEEHHQSLDPNCPRDLTDCLLVEMEKEKHSAERLYTMDGITVTVAKLFFAGTETTSTTLRYGLLILMKYPEIEEKDHEEIDRPIGPSRIPAIKDRQEMPYMDAVVHEIQRFITLVPSNCPHEATRQTIFRGYLIPKETVVVPTLDSVLYDNQEFPDPEKFKPEEFLNENGKFKYSDYFKPFSTGKRVCAGEGLARMELFLLLCAILQHFNLKPLVDPKDIDLSPRHIGHGCIPPRYKLEVIPRS*"

dna = Seq(mutated, generic_dna)
gen_0 = str(dna.translate())

dna2 = Seq(original, generic_dna) #dna2 :)
gen_2 = str(dna2.translate())

#inmutable sites, warning: one-based numbering
#estaticos = set(range(233,244)+[105,210,295,364,368,478]+range(1,22)+range(302,306)+[430,437,446])
estaticos = set(range(233,244)+range(1,22)+range(302,306)+[430,437,446])

# CXPB  is the probability with which two individuals
#       are crossed
#
# MUTPB is the probability for mutating an individual
#
# NGEN  is the number of generations for which the
#       evolution runs
# PSIZE is population size

CXPB, MUTPB, NGEN, PSIZE = 1.0, 0.9, 500, 200

#genamac = partial(random.choice, "ARNDCEQGHILKMFPSTWYV")

blosum62 ={
"A": { "R": -1, "N": -2, "D": -2, "C":  0, "Q": -1, "E": -1, "G":  0, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S":  1, "T":  0, "W": -3, "Y": -2, "V":  0 },
"R": { "A":-1, "N":  0, "D": -2, "C": -3, "Q":  1, "E":  0, "G": -2, "H":  0, "I": -3, "L": -2, "K":  2, "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3 },
"N": { "A":-2, "R":  0, "D":  1, "C": -3, "Q":  0, "E":  0, "G":  0, "H":  1, "I": -3, "L": -3, "K":  0, "M": -2, "F": -3, "P": -2, "S":  1, "T":  0, "W": -4, "Y": -2, "V": -3 },
"D": { "A":-2, "R": -2, "N":  1, "C": -3, "Q":  0, "E":  2, "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3, "F": -3, "P": -1, "S":  0, "T": -1, "W": -4, "Y": -3, "V": -3 },
"C": { "A": 0, "R": -3, "N": -3, "D": -3, "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2, "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1 },
"Q": { "A":-1, "R":  1, "N":  0, "D":  0, "C": -3, "E":  2, "G": -2, "H":  0, "I": -3, "L": -2, "K":  1, "M":  0, "F": -3, "P": -1, "S":  0, "T": -1, "W": -2, "Y": -1, "V": -2 },
"E": { "A":-1, "R":  0, "N":  0, "D":  2, "C": -4, "Q":  2, "G": -2, "H":  0, "I": -3, "L": -3, "K":  1, "M": -2, "F": -3, "P": -1, "S":  0, "T": -1, "W": -3, "Y": -2, "V": -2 },
"G": { "A": 0, "R": -2, "N":  0, "D": -1, "C": -3, "Q": -2, "E": -2, "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S":  0, "T": -2, "W": -2, "Y": -3, "V": -3 },
"H": { "A":-2, "R":  0, "N":  1, "D": -1, "C": -3, "Q":  0, "E":  0, "G": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y":  2, "V": -3 },
"I": { "A":-1, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -3, "E": -3, "G": -4, "H": -3, "L":  2, "K": -3, "M":  1, "F":  0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V":  3 },
"L": { "A":-1, "R": -2, "N": -3, "D": -4, "C": -1, "Q": -2, "E": -3, "G": -4, "H": -3, "I":  2, "K": -2, "M":  2, "F":  0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V":  1 },
"K": { "A":-1, "R":  2, "N":  0, "D": -1, "C": -3, "Q":  1, "E":  1, "G": -2, "H": -1, "I": -3, "L": -2, "M": -1, "F": -3, "P": -1, "S":  0, "T": -1, "W": -3, "Y": -2, "V": -2 },
"M": { "A":-1, "R": -1, "N": -2, "D": -3, "C": -1, "Q":  0, "E": -2, "G": -3, "H": -2, "I":  1, "L":  2, "K": -1, "F":  0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V":  1 },
"F": { "A":-2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I":  0, "L":  0, "K": -3, "M":  0, "P": -4, "S": -2, "T": -2, "W":  1, "Y":  3, "V": -1 },
"P": { "A":-1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2 },
"S": { "A": 1, "R": -1, "N":  1, "D":  0, "C": -1, "Q":  0, "E":  0, "G":  0, "H": -1, "I": -2, "L": -2, "K":  0, "M": -1, "F": -2, "P": -1, "T":  1, "W": -3, "Y": -2, "V": -2 },
"T": { "A": 0, "R": -1, "N":  0, "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S":  1, "W": -2, "Y": -2, "V":  0 },
"W": { "A":-3, "R": -3, "N": -4, "D": -4, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F":  1, "P": -4, "S": -3, "T": -2, "Y":  2, "V": -3 },
"Y": { "A":-2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H":  2, "I": -1, "L": -1, "K": -2, "M": -1, "F":  3, "P": -3, "S": -2, "T": -2, "W":  2, "V": -1 },
"V": { "A": 0, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I":  3, "L":  1, "K": -2, "M":  1, "F": -1, "P": -2, "S": -2, "T":  0, "W": -3, "Y": -1 },
"*": { "*": 0 }
}

cache = {}

def evalOneMax(individual):
   seq = "".join(individual)
   if seq in cache:
      return cache[seq]
   #5% from the original sequence, max
   elif sum(a!=b for a,b in zip(gen_2, seq)) > 0.05 * len(gen_2):
      return 0, 0
   else:
      #check solubility
      print seq
      while True:
         try:
            r = requests.post("http://cblab.meiyaku.jp/ESPRESSO/CheckInput.php", 
                              data={'System': 'Ecoli', 'Model': 'SOL', 'Data': 'AA', 'Sequence': seq[:-1]}, timeout=10)
            calcID = re.search('value="([^"]+)" name="calcID"',r.text).group(1)
            r = requests.post("http://cblab.meiyaku.jp/ESPRESSO/Prediction.php", 
                              data={'System': 'Ecoli', 'Model': 'SOL', 'Data': 'AA', 'Sequence': seq[:-1], 'calcID': calcID}, timeout=10)
            break
         except:
            print "Oops! Try again..."
      

      url = "http://cblab.meiyaku.jp/ESPRESSO/RESULTS/" + str(calcID).strip() + ".html"
      print url
      time.sleep(40)
      c=0
      while c < 60*3:
         while True:
            try:
               r = requests.post(url, timeout=10)
               break
            except:
               print "Oops! Try again..."
         if re.search("Insoluble \((.+)\)", r.text) or re.search("Soluble \((.+)\)", r.text):
           break
         time.sleep(1)
         c += 1
      if re.search("Insoluble \((.+)\)", r.text):
         solubility = 1 - float(re.search("Insoluble \((.+)\)", r.text).group(1))
      elif re.search("Soluble \((.+)\)", r.text):
         solubility = float(re.search("Soluble \((.+)\)", r.text).group(1))
      else:
         solubility = 0
         print "ERROR"
      #Optimizer, to find DNA sequence
      while True:
         try:
            r = requests.post("http://genomes.urv.es/OPTIMIZER/obtimized.php", data = {"PROTop": "Robust", "invert": "",
                              "sequence": seq[:-1], "type": "Proteins", "weigthOP": "optimizer", "format": "CodonUsageDB",
                              "typeRef": "codon", "codonWi": "", "nom[0]": "ecoli", "geneticcode": "11"}, timeout=10)
            dnaseq = re.search("name=opsequence value=([ACGT]+)>", r.text).group(1)
            break
         except:
            print "Oops! Try again..."
      #check expression
      while True:
         try:
            r = requests.post("http://cblab.meiyaku.jp/ESPRESSO/CheckInput.php", 
                              data={'System': 'Ecoli', 'Model': 'EXP', 'Data': 'nt', 'Sequence': dnaseq}, timeout=10)
            calcID = re.search('value="([^"]+)" name="calcID"',r.text).group(1)
            r = requests.post("http://cblab.meiyaku.jp/ESPRESSO/Prediction.php", 
                              data={'System': 'Ecoli', 'Model': 'EXP', 'Data': 'nt', 'Sequence': dnaseq, 'calcID': calcID}, timeout=10)
            break
         except:
            print "Oops! Try again..."
      url = "http://cblab.meiyaku.jp/ESPRESSO/RESULTS/" + str(calcID).strip() + ".html"
      print url
      time.sleep(40)
      c=0
      while c < 60*3:
         while True:
            try:
               r = requests.post(url)
               break
            except:
               print "Oops! Try again..."
         if re.search(">Not Expressed \((.+)\)", r.text) or re.search(">Expressed \((.+)\)", r.text):
           break
         if re.search("Your sequence cannot be translated into protein sequence.", r.text):
           break
         time.sleep(1)
         c += 1
      if re.search(">Not Expressed \((.+)\)", r.text):
         expression = 1 - float(re.search(">Not Expressed \((.+)\)", r.text).group(1))
      elif re.search(">Expressed \((.+)\)", r.text):
         expression = float(re.search(">Expressed \((.+)\)", r.text).group(1))
      elif re.search("Your sequence cannot be translated into protein sequence.", r.text):
         expression = 0
      else:
         expression = 0
         print "ERROR"
   cache[seq] = (solubility, expression)
   print cache[seq]
   return cache[seq]
   
#def evalOneMax(individual):
#   return individual.count("V"),

def my_mutation(individual, indpb=0):
   while True:
      i = random.randint(0,len(individual)-1)
      if i+1 not in estaticos:
         break
   options_mut = list(chain(*[repeat(k,x - min(blosum62[individual[i]].values()) + 1) for k,x in blosum62["A"].items()]))
   individual[i] = random.choice(options_mut)
   return individual,

def my_cross(dad, mom):
   i = random.randint(1,len(dad)-2)
   son = creator.Individual(dad[:i] + mom[i:])
   daughter = creator.Individual(mom[:i] + dad[i:])
   return son, daughter

ngen_0 = partial(next, cycle(gen_0))
ngen_2 = partial(next, cycle(gen_2))
def gen_population():
   #original
   if random.random() < 1.00:
     ngen_2()
     return ngen_0()
   #mutada
   ngen_0()
   return ngen_2() 

creator.create("FitnessMulti", base.Fitness, weights=(1.0, 1.0))
creator.create("Individual", list, fitness=creator.FitnessMulti)

toolbox = base.Toolbox()
toolbox.register("attr_gen", gen_population)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_gen, len(gen_0))
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", evalOneMax)
toolbox.register("mate", my_cross)
toolbox.register("mutate", my_mutation, indpb=0.5)
toolbox.register("select", tools.selTournament, tournsize=2)

pop = toolbox.population(n=PSIZE)
    
print("Start of evolution")
    
# Evaluate the entire population
fitnesses = list(map(toolbox.evaluate, pop))
for ind, fit in zip(pop, fitnesses):
   ind.fitness.values = fit
    
print("  Evaluated %i individuals" % len(pop))
    
# Begin the evolution
for g in range(NGEN):
   print("-- Generation %i --" % g)

   # Select the next generation individuals
   offspring = toolbox.select(pop, len(pop))
   # Clone the selected individuals
   offspring = list(map(toolbox.clone, offspring))
    
   # Apply crossover and mutation on the offspring
   for child1, child2 in zip(offspring[::2], offspring[1::2]):
      # cross two individuals with probability CXPB
      if random.random() < CXPB:
         toolbox.mate(child1, child2)
         # fitness values of the children
         # must be recalculated later
         del child1.fitness.values
         del child2.fitness.values

   for mutant in offspring:
      # mutate an individual with probability MUTPB
      if random.random() < MUTPB:
         toolbox.mutate(mutant)
         del mutant.fitness.values
    
   # Evaluate the individuals with an invalid fitness
   invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
   fitnesses = map(toolbox.evaluate, invalid_ind)
   for ind, fit in zip(invalid_ind, fitnesses):
      ind.fitness.values = fit
        
   print("  Evaluated %i individuals" % len(invalid_ind))
        
   # The population is entirely replaced by the offspring
   pop[:] = offspring
        
   # Gather all the fitnesses in one list and print the stats
   fits = [ind.fitness.values[0] for ind in pop]
        
   length = len(pop)
   mean = sum(fits) / length
   sum2 = sum(x*x for x in fits)
   std = abs(sum2 / length - mean**2)**0.5
   
   print("".join(tools.selBest(pop, 1)[0]))     
   print("  Min %s" % min(fits))
   print("  Max %s" % max(fits))
   print("  Avg %s" % mean)
   print("  Std %s" % std)
    
best_ind = tools.selBest(pop, 1)[0]
print("Best individual is %s, %s" % ("".join(best_ind), best_ind.fitness.values))

