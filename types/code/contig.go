// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"

	"github.com/biogo/examples/contig"
)

type offsetSeq struct {
	seq    seq.Sequence
	offset int
}

func main() {
	seqs := []offsetSeq{
		{
			seq:    linear.NewSeq("id1", alphabet.BytesToLetters([]byte("AGTC")), alphabet.DNA),
			offset: 2,
		},
		{
			seq:    linear.NewSeq("id2", alphabet.BytesToLetters([]byte("ACGT")), alphabet.DNA),
			offset: 15,
		},
	}

	con, err := contig.New("super contig", 1, alphabet.DNA)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	con.Relaxed(true)

	for _, os := range seqs {
		os.seq.SetOffset(os.offset)
		con.Insert(os.seq)
	}

	wrapLength := 10

	fmt.Printf("%*a\n", wrapLength, con)
	con.RevComp()
	fmt.Printf("%*a\n", wrapLength, con)
}
