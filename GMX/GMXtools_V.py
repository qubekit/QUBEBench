import GMXmdp


GMXmdp.create_top('UNK')
GMXmdp.setup_BOX('UNK')
GMXmdp.minimize_steep(resname='UNK')
GMXmdp.NVT_Equilibrate(resname='UNK')
GMXmdp.NPT_Equilibrate(resname='UNK')
GMXmdp.NPT_Production(resname='UNK')

# for i in range(1):
#     GMXmdp.minimize_steep(resname='UNK', window=i, anh='Q')
#     GMXmdp.NVT_Equilibrate(resname='UNK', window=i, anh='Q', solname='SOL')
#     GMXmdp.NPT_Equilibrate(resname='UNK', window=i, anh='Q', solname='SOL')
#     GMXmdp.NPT_Production(resname='UNK', window=i, anh='Q', solname='SOL')
#     GMXmdp.write_fep_job(window=i)

for i in range(16):
    GMXmdp.minimize_steep(resname='UNK', window=i, anh='Vdw')
    GMXmdp.NVT_Equilibrate(resname='UNK', window=i, anh='Vdw', solname='Water')
    GMXmdp.NPT_Equilibrate(resname='UNK', window=i, anh='Vdw', solname='Water')
    GMXmdp.NPT_Production(resname='UNK', window=i, anh='Vdw', solname='Water')
    GMXmdp.write_fep_job(window=i)
