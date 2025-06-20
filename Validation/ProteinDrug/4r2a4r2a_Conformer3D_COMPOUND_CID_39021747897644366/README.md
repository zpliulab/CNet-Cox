Output files of CB-Dock2 server:

For "Template-based Blind Docking", if some templates were found in our database, you will get some output files that start with "fitdock".
	
	1. Docking conformation of the query ligand by FitDock: 

		fitdock_em_<num>.<score>.mol2 
			- MOL2 format of the query ligand conformation after docking
		fitdock_em_<num>.<score>.pdb
			- PDB format of the query ligand conformation after docking
		fitdock_em_<num>.<score>.complex.pdb
			- PDB format of the complex conformation after docking

		Nomenclature: the '<num>' refers to the template id and the '<score>' refers to the affinity score of the docking conformation.

	2. Conformation of template ligands 

		fitdock_lt_<num>.mol2
			- MOL2 format of the template ligand 

		Nomenclature: the '<num>' refers to the template id

	3. Conformation of template proteins 

		fitdock_pt_<num>.pdb
			- PDB format of the template protein 

		Nomenclature: the '<num>' refers to the template id

For "Structure-based Blind Docking", you will get some output files that start with the name of query protein.

	4. Docking conformation of the query ligand by AutoDock Vina:

		<protein_name>:<ligand_name>_out_<num>.<score>.mol2
			- MOL2 format of the query ligand conformation after docking
		<protein_name>:<ligand_name>_out_<num>.<score>.pdb
			- PDB format of the query ligand conformation after docking
		<protein_name>:<ligand_name>_out_<num>.<score>.complex.pdb
			- PDB format of the complex conformation after docking

		Nomenclature: the '<protein_name>' refers to the name of query protein, the '<ligand_name>' refers to the name of query ligand, the '<num>' refers to the CurPocket id and the '<score>' refers to the affinity score of the docking conformation.
	
