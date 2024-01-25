python Recortar.py
chmod 777 dopsfex_mask.sh
./dopsfex_mask.sh
python inputs_galfitm.py
chmod 777 galfit_215_STRIPE82-0073.input
./galfitm-1.4.4-linux-x86_64 inputs/galfit_215_STRIPE82-0073.input
python leer_outputs.py
