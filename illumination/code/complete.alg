Function Collision(R reads):
	Let S be a spatially indexed store of reads
	Let C and M be 0
	For each read pair r in R:
		Increment M
		Store r in S
	For each read pair r in R:
		Let c be closest read to r in S
		If mapping-of(c) overlaps mapping-of(r)
			Increment C
	Return C, M.
