
# Figure 12 
./vipss -i ../data/hand_ok/input.xyz -l 0 -s 200

# Figure 1 
./vipss -i ../data/walrus/input.xyz -l 0.003 -s 100

# Figure 9, 11 
./vipss -i ../data/bathtub/input.xyz -l 0 -s 200

# Figure 12 
./vipss -i ../data/phone/input.xyz -l 0 -s 100

# Figure 5
./vipss -i ../data/planck/multisample_n500/input.xyz -l 0 -s 100 -t
./vipss -i ../data/planck/multisample_n1000/input.xyz -l 0 -s 100 -t
./vipss -i ../data/planck/multisample_n2000/input.xyz -l 0 -s 100 -t
./vipss -i ../data/planck/multisample_n4000/input.xyz -l 0 -s 100 -t

# Figure 6
./vipss -i ../data/wireframes/doghead/input.xyz -l 0 -s 100
./vipss -i ../data/wireframes/phone/input.xyz -l 0 -s 100
./vipss -i ../data/wireframes/trebol/input.xyz -l 0 -s 100
./vipss -i ../data/wireframes/spring/input.xyz -l 0 -s 100
./vipss -i ../data/surfaces_500/vertebra/input.xyz -l 0 -s 100
./vipss -i ../data/surfaces_500/hand/input.xyz -l 0 -s 100

# Figure 4
./vipss -i ../data/torus/crosscut/input.xyz -l 0 -s 50
./vipss -i ../data/torus/halfsamplel500_r25/input.xyz -l 0 -s 50
./vipss -i ../data/torus/multisample_n25/input.xyz -l 0 -s 50
./vipss -i ../data/torus/multisample_n50/input.xyz -l 0 -s 50
./vipss -i ../data/torus/wires/input.xyz -l 0 -s 50
./vipss -i ../data/torus/parallelcut/input.xyz -l 0 -s 50

# Figure 7
./vipss -i ../data/noise_kitten/kitten_h0.01/input.xyz -o ../data/noise_kitten/kitten_h001/l0001_ -l 0.001 -s 100
./vipss -i ../data/noise_kitten/kitten_h0.01/input.xyz -o ../data/noise_kitten/kitten_h001/l001_ -l 0.01 -s 100

./vipss -i ../data/noise_kitten/kitten_h0.04/input.xyz -o ../data/noise_kitten/kitten_h004/l0001_ -l 0.001 -s 100
./vipss -i ../data/noise_kitten/kitten_h0.04/input.xyz -o ../data/noise_kitten/kitten_h004/l001_ -l 0.01 -s 100
