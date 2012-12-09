Function PILER-DF-NEW(P piles, t integer):
	Let G be a graph where nodes represent repeat instances
	For each pile p in P:
		Let C = clusterise(p)
		For each cluster c in C:
			For each image Q in c:
				Let T = partner(Q)
				Add edge c—cluster-of(T) to G // HL
	Return connected components of G of order ≥ t.