Function Collision(R reads):
	Let S be a spatially indexed store of reads
	Let C and D be 0
	For each read pair r in R:
		Increment D
		If is-discordant(r):
			Store r in S
	For each read pair r in R:
		If is-not-discordant(r):
			Let c be closest read to r in S
			If mapping-of(c) overlaps mapping-of(r)
				Increment C
	Return C, D.
