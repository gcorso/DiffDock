# Makefile for aspects of installing DiffDock

# Pip install and copy files
# May need root or bmaps-server permission
install: DEST_DIR=/opt/diffdock
install:
	-mkdir $(DEST_DIR)
	cd $(DEST_DIR) && python3 -m venv venv
	rsync -az --exclude=.git --exclude=coniferpoint ../* $(DEST_DIR)/DiffDockCP
	# Make sure we fail with an error before invoking the foreach. If venv fails, the foreach will still attempt all the installs.
	cd $(DEST_DIR) && . venv/bin/activate
	# Use foreach to install each python package individually because it fails when it tries to
	# do them all together.
	cd $(DEST_DIR) && . venv/bin/activate && $(foreach pkg, $(shell cat ./requirements.txt), pip install $(pkg);)

# Download and cache a large torch model
seed: DIFFDOCK_DIR=/opt/diffdock
seed:
	cd $(DIFFDOCK_DIR)/DiffDockCP && . ../venv/bin/activate && TORCH_HOME=$(DIFFDOCK_DIR) python -m inference

# Test run
run: DIFFDOCK_DIR=/opt/diffdock
run: OUT_DIR=
run: PDB=
run: SDF=
run:
	if test -z $(OUT_DIR) || test -z $(PDB) || test-z $(SDF); then echo "Please define PDB SDF and OUT_DIR"; exit 1; fi
	cd $(DIFFDOCK_DIR)/DiffDockCP && . ../venv/bin/activate && TORCH_HOME=$(DIFFDOCK_DIR) python -m inference --complex_name 'TEST' --protein_path $(PDB) --ligand_description $(SDF) --out_dir $(OUT_DIR)

# Create the root directory (may need root permission)
root: DEST_DIR=/opt/diffdock
root: OWNER=bmaps-server:bmaps-server
root:
	-mkdir -p $(DEST_DIR)
	chown -R $(OWNER) $(DEST_DIR)
