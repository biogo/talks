// Copyright ©2012 Dan Kortschak <dan.kortschak@adelaide.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

package main

import (
	"bytes"
	"code.google.com/p/biogo.external/muscle"
	"code.google.com/p/biogo/exp/alphabet"
	"code.google.com/p/biogo/exp/seq"
	"code.google.com/p/biogo/exp/seq/linear"
	"code.google.com/p/biogo/exp/seq/multi"
	"code.google.com/p/biogo/exp/seqio/fasta"
	"fmt"
	"io"
	"strings"
)

var s = `>71.2259 lcl|scaffold_41:8288143+
CCCCAAATTCTCATAAAAAGACCAGACTTAATGGTCTGACTGAGACTAGAGGAATCCCGG
TGGTCATGGTCCCCAAACCTTCTGTTGGCCCAGGACAGGAACCATTCCCGAAGACAACTC
ATCAGACACGGAAGGGACTGGACAATGGGTAGGAGAGAGATGCTGACGAAGAGTGAGCTA
CTTGTATCAGGTGGACACTTGAGACTGTGTTGGCATCTCCTGTCTGGAGGGGAGATAGGA
GGGTAGAGAGGGTTAGAAACTGGCAAAATCGTCATGAAAGGAGGGACTGGAAGGAGGGAG
CGGGCTGACTCAGTAGGGGGAGAGTAAGTGGGAGTATGGAGTAAGGTGTATATAAGCTTA
TATGTGACAGATTGACTTGATTTGTAAACTTTCACTTAAAGCACAATAAAAATTATTTTT
TAAAAAATTGTTT
>71.2259 lcl|scaffold_41:11597466-
ATTATTATTTTTTTAAATAATTTTTATTGTGTTTTAAGGGAAAGTTTGCAAATCAAGTCA
GTCTCTCACATATAACCTTATATACACCTTACTCCATACTCCCATTTACTCTCCCCCTAA
TGAGTCAGCCCGCTCCCTCCTTCCGGTCTCTCCTTTCTTGACGATTTTGTCAGTTTCTAA
CCCTCTCTACCCTTCTATCTCTCCTCCAGACAGGAGATGCCAACACTGTCTCAAGTGTCC
ACTTGATACAAGTAGCTCACTCTTCGTCAGCATCTCTCTCCAACCCATTGTCCAGTCCCT
GCCATGTCTGATGAGTTGTCTTTGGGAATGGTTCCTGTCCTGGGCCAACAGAAGGTTTGG
GGACCATGACCGCTGGGATTCCTCTAGTCTCAGTCAGACCATTAAGTCTGGTCTTTTTAT
GAGA
>71.2259 lcl|scaffold_45:2724255+
ATAAAAAGACCAGACTTAATGGTCTGACTGAGACTAGAAGAATCCCGGTGGCCATGGTCC
CCAAACCTTCTGTTGGCCCAGGACAGGAACCATTCCCGAAGACAATTCATCAGACATGGA
AGGGACTGGACAATGGGTTGGAGAGAGATGCTGATAAAGAGTGAGCTACTTGTATCAGGT
GGACGTTTGAGACTGTATTGGCATCTCCTGTCTGGAGGGGAGATAGGGTAGAGAGGGTTA
GAAACTGGCAAAACGGTCACGAAAGGAGAGACTGGAAGAAGGGAGCAGGCTGACTCATTA
GGGGGAGAGTAAATGGGAGTATGTAGTAAGGTGTATATAAGCTTACATGTGACAGACTGA
CTTGATTTGTAAACTTTCACTTAAAGCACAATAAAAATTATTTTTTAAAAATTTGCC
`

func main() {
	m, err := muscle.Muscle{Quiet: true}.BuildCommand()
	if err != nil {
		panic(err)
	}
	m.Stdin = strings.NewReader(s)
	m.Stdout = &bytes.Buffer{}
	m.Run()
	// {CONS OMIT
	var (
		r = fasta.NewReader(m.Stdout.(io.Reader), &linear.Seq{
			Annotation: seq.Annotation{Alpha: alphabet.DNA},
		})
		ms = &multi.Multi{
			ColumnConsense: seq.DefaultQConsensus,
		}
	)
	for {
		s, err := r.Read()
		if err != nil {
			break
		}
		ms.SetName(s.Name())
		ms.Add(s)
	}
	c := ms.Consensus(false)
	c.Threshold = 42
	c.QFilter = seq.CaseFilter
	fmt.Printf("%60a\n", c)
	// CONS} OMIT
}
