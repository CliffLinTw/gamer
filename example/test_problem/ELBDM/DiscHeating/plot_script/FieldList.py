import yt
ds = yt.load("/work1/clifflin/gamer-fork/bin/Plummer/Data_000000")
for i in sorted(ds.field_list):
	print(i)
