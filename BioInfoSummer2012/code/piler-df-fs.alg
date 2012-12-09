Function PILER-DF(P piles, t integer):
	Let G be a graph where nodes represent repeat instances
	For each pile p in P:
		For each image Q in p:
			Let T = partner(Q)
			If is-global-image(Q) and is-global-image(T):
				Add edge p—pile-of(T) to G
	Return connected components of G of order ≥ t. // HL