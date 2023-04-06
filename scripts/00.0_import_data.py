import os

CHROMOPLASTIC_PROJECT_PATH = "/mnt/etemp/luca/Projects/WorkInProgress/ChromoPlastiC"
COMPARTMENTS_PATH = "data/calder"


os.makedirs(COMPARTMENTS_PATH, exist_ok=True)

print("Importing Calder calls from ChromoPlastiC project")
CHROMOPLASTIC_HIC_PATH = os.path.join(CHROMOPLASTIC_PROJECT_PATH, "data", "processed_data", "hic")
for sample in os.listdir(CHROMOPLASTIC_HIC_PATH):
	sample_calder_path = os.path.join(CHROMOPLASTIC_HIC_PATH,
									  sample,
									  "compartments",
									  "calder",
									  "sub_compartments",
									  "all_sub_compartments.bed")
	sample_out_path = os.path.join(COMPARTMENTS_PATH, f"{sample}_CompartmentDomains.bed")
	if os.path.isfile(sample_calder_path) and (not os.path.exists(sample_out_path)):
		print(f"- {sample}")
		os.symlink(sample_calder_path, sample_out_path)
