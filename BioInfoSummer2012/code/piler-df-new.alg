Function PILER-DF-NEW(P piles, t integer):
	Let G be a graph where nodes represent repeat instances
	For each pile p in P:
		Let C = clusterise(p) // HL
		For each cluster c in C: // HL
			For each image Q in c: // HL
				Let T = partner(Q)
				If is-global-image(Q) and is-global-image(T):
					Add edge c—cluster-of(T) to G // by analogy to pile-of() // HL
	Return connected components of G of order ≥ t.