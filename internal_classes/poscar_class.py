import numpy as np
import os

def vesta_to_coord(vesta_string):
    return [np.array(vesta_string.split(), dtype=float)]

class poscar:
    def __init__(self, location="./POSCAR"):
        with open(location) as file:
            self.system = file.readline().strip()
            self.factor = file.readline().strip()
            self.a = file.readline().strip()
            self.b = file.readline().strip()
            self.c = file.readline().strip()
            self.species = file.readline().split()
            self.species_count = np.fromstring(file.readline().strip(), dtype=int, sep=' ')
            
            self.seldyn = file.readline().strip()
            if self.seldyn[0] in ['s', 'S']:
                self.direct = file.readline().strip()
        
            if self.seldyn[0] not in ['s', 'S']:
                print("Adding seldyn line")
                self.direct = self.seldyn
                self.seldyn = 'seldyn_appended_by_script'
            self.total_atoms = self.species_count.sum()
            self.coordinates = []
            for i in range(self.total_atoms):
                self.coordinates.append(file.readline().strip())
            self.coordinates_numbers = [coordinate.split()[0:3] for coordinate in self.coordinates]
            self.coordinates_numbers = np.array(self.coordinates_numbers, dtype=float)
            
            self.coordinates_seldyn = [' '.join(coordinate.split()[3:6]) for coordinate in self.coordinates]
            if all([len(x.split()) == 3 for x in self.coordinates_seldyn]):
                self.coordinates_seldyn = np.array(self.coordinates_seldyn)
            else:
                print("Adding True values for seldyn to all coordinates")
                self.coordinates_seldyn = []
                for i in range(self.total_atoms):
                    self.coordinates_seldyn.append("T T T") 
                self.coordinates_seldyn = np.array(self.coordinates_seldyn)

    def save(self, location="./POSCAR_mod"):
        if os.path.exists(location):
            print("Overwriting file")
            os.unlink(location)
        with open(location, 'w') as file:
            np.set_printoptions(suppress=True)
            file.write(self.system+'\n')
            file.write(self.factor+'\n')
            file.write(self.a+'\n')
            file.write(self.b+'\n')
            file.write(self.c+'\n')
            file.write(" ".join(self.species)+'\n')
            file.write(' '.join(map(str, self.species_count))+'\n')
            file.write(self.seldyn+'\n')
            file.write(self.direct+'\n')
            for i in range(self.total_atoms):
                file.write(' '.join(map(str, self.coordinates_numbers[i])) + ' ' + self.coordinates_seldyn[i] + '\n')
                
    def unfreeze(self):
        for i in range(self.total_atoms):
             self.coordinates_seldyn[i] = 'T T T'
                
    def freeze(self, threshold):
        for i in range(self.total_atoms):
            if self.coordinates_numbers[i][2] < threshold:
                self.coordinates_seldyn[i] = 'F F F'
            else:
                self.coordinates_seldyn[i] = 'T T T'
    
    def add_species(self, species, count, coord):
        if len(coord) != count:
            print("Warning! Number of coordinates does not correspond to number of added species")
        if self.species[-1] == species:
            self.species_count[-1] = self.species_count[-1] + count
        else:
            self.species = np.append(self.species, species)
            self.species_count = np.append(self.species_count, count)

        self.coordinates_numbers = np.append(self.coordinates_numbers, coord, axis = 0)
        self.coordinates_seldyn = np.append(self.coordinates_seldyn, np.repeat(['T T T'], count, axis = 0), axis = 0)
        self.total_atoms = self.total_atoms + count

                
    
                           