# Conifer Point DiffDock Setup

(intent: invoke `make root` as root to create /opt/diffdock, then invoke `make install` as bmaps-server user to make the child directory)

1. Identify an owner user (or user:group) and a destination folder
2. `sudo make root OWNER=<OWNER USER|USER:GROUP>` 
3. `sudo -u <OWNER USER> make install DEST_DIR=<DESTINATION FOLDER>`
4. `sudo -u <OWNER USER> make seed DIFFDOCK_DIR=<DESTINATION FOLDER>`
