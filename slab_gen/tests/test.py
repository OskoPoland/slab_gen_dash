import json
from slab_gen import InterfaceBuilder, Slab

if __name__ == "__main__":
    # -- Creating a new slab
    test_slab_path = 'C:\\Users\Vinr5\\OneDrive\\Documents\\Repositories\\Adelstein-research\\codebase\\slab_gen-Arye\\scripts\\Build_interface\\Structures\\Li2CO3.POSCAR.vasp'
    test_slab_path2 = 'C:\\Users\\Vinr5\\OneDrive\\Documents\\Repositories\\Adelstein-research\\codebase\\slab_gen-Arye\\scripts\\Build_interface\\Structures\\LiF.POSCAR.vasp'
    test_out_path = 'C:\\Users\\Vinr5\\OneDrive\\Documents\\Repositories\\Adelstein-research\\codebase\\slab_gen-Arye\\slab_gen\\tests'
    slab1 = Slab(test_slab_path, (0, 0, 1)) # nlay and vacheight are defaults 
    slab2 = Slab(test_slab_path2, (0, 0, 1))

    # -- Slab serialization for easy reuse
    slab1.to_json(name="slab1_test", out=test_out_path)
    slab1_copy = Slab.from_json(path=test_out_path + '\\slab1_test.json')
    slab1_copy.to_json(name='slab1_test_copy', out=test_out_path)
   
    # -- Generating an interface
    #       - Interface (only 1 - not an optimal one) is fully generated when the constructor is called
    #       - Previous gen_candidates parameters are now passed into the constructor
    #       - Every parameters has a getter and setter for easier parameter modification after generation
    test_interface = InterfaceBuilder(slab1, slab2, ncand=10, nmax=6, mmax=10)
        
    # -- to generate translations
    #       - Will populate the self.interface list with all generated candidates and sort them by strain
    #       - List of interfaces is emptied before each generation 
    #       - If you use the getters and setters to modify params, simply calling generate translations
    #         will overwrite previous generated translations
    test_interface.generate_translations()
    
    # -- to write the translation .dat files specifically mesh, prim_cell, and super_cell
    test_interface.write_translations(output=test_out_path)
    
    # -- to write all of the candidates
    test_interface.write_candidates(output=test_out_path)
    
    # -- to write all the interfaces
    test_interface.write_interfaces(output=test_out_path)
    
    
    
    
    

    
