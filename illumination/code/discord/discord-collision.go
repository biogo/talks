// Copyright ©2013 The bíogo.talks Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/biogo/boom"
	"github.com/biogo/illumina"
	"github.com/biogo/store/kdtree"
	"fmt"
	"math"
	"os"
)

// Offset definition for overlap comparison.
type offset struct {
	dist  int
	label string
}

// The offsets we are interested in.
var offsets = []offset{
	{0, "Coincide"},
	{1e2, "Adjacent"},
	{1e3, "At1k"},
	{1e4, "At10k"},
}

// This allows analysis of systems where there is some obliquity.
const (
	xunit = 37.5 // The width of a coordinate in nm.
	yunit = 37.5 // The height of a coordinate in nm.
)

// randoms is the number of random points to choose when determining median points
// for pivot operations.
var randoms = 100

// boomIllumina wraps boom.Record in order to satisfy illumina.Interface.
type boomIllumina struct{ *boom.Record }

func (b boomIllumina) Description() string { return "" }

// mapping is a terse representation of bam mapping data.
type mapping struct {
	Segment    string
	Start, End int
}

// tileAddress is a hashable unique tile identifier.
type tileAddress struct {
	FlowCell string
	Lane     int8
	Tile     int
}

// illuminaRecord stores mapping an illumina meta data and satisfies the kdtree.Comparable
// interface using the flow cell coordinates as the point data.
type illuminaRecord struct {
	A, B mapping
	illumina.Metadata
}

// illuminaRecord} OMIT

// newRecord returns an illumina record based on two boom.Records and a set of reference names.
func newRecord(r [2]*boom.Record, names []string) (*illuminaRecord, error) {
	m, err := illumina.Parse(boomIllumina{r[0]}) // They are a pair, so we only parse one.
	if err != nil {
		return nil, err
	}
	return &illuminaRecord{
		A: mapping{
			Segment: names[r[0].RefID()],
			Start:   r[0].Start(),
			End:     r[0].End(),
		},
		B: mapping{
			Segment: names[r[1].RefID()],
			Start:   r[1].Start(),
			End:     r[1].End(),
		},
		Metadata: m,
	}, nil
}

// overlap returns true if there is any overlap between the reads of the two provided
// illuminaRecords, after applying the at offset to a.
func overlap(a, b *illuminaRecord, at int) bool {
	return (a.A.Segment == b.A.Segment && a.A.End-at > b.A.Start && a.A.Start-at < b.A.End) ||
		(a.B.Segment == b.B.Segment && a.B.End-at > b.B.Start && a.B.Start-at < b.B.End) ||
		(a.A.Segment == b.B.Segment && a.A.End-at > b.B.Start && a.A.Start-at < b.B.End) ||
		(a.B.Segment == b.A.Segment && a.B.End-at > b.A.Start && a.B.Start-at < b.A.End)
}

// satisfy the kdtree.Comparable interface. The dimensions are:
//  0 = X
//  1 = Y
// {illuminaRecord methods OMIT
func (p *illuminaRecord) Compare(c kdtree.Comparable, d kdtree.Dim) float64 {
	q := c.(*illuminaRecord)
	switch d {
	case 0:
		return float64(p.Coordinate.X-q.Coordinate.X) * xunit
	case 1:
		return float64(p.Coordinate.Y-q.Coordinate.Y) * yunit
	default:
		panic("illegal dimension")
	}
}
func (p *illuminaRecord) Dims() int { return 2 }
func (p *illuminaRecord) Distance(c kdtree.Comparable) float64 {
	q := c.(*illuminaRecord)
	x := float64(p.Coordinate.X-q.Coordinate.X) * xunit
	y := float64(p.Coordinate.Y-q.Coordinate.Y) * yunit
	return x*x + y*y
}

// illuminaRecord methods} OMIT

// illuminaRecords is a collection of the illuminaRecord type that satisfies kdtree.Interface.
// This type is implemented to improve performance of queries in the second pass.
type illuminaRecords []*illuminaRecord

func (p illuminaRecords) Index(i int) kdtree.Comparable { return p[i] }
func (p illuminaRecords) Len() int                      { return len(p) }
func (p illuminaRecords) Pivot(d kdtree.Dim) int {
	return plane{illuminaRecords: p, Dim: d}.Pivot()
}
func (p illuminaRecords) Slice(start, end int) kdtree.Interface { return p[start:end] }

// illuminaRecords} OMIT

// plane is required to help illuminaRecords.
type plane struct {
	kdtree.Dim
	illuminaRecords
}

func (p plane) Less(i, j int) bool {
	switch p.Dim {
	case 0:
		return p.illuminaRecords[i].Coordinate.X < p.illuminaRecords[j].Coordinate.X
	case 1:
		return p.illuminaRecords[i].Coordinate.Y < p.illuminaRecords[j].Coordinate.Y
	default:
		panic("illegal dimension")
	}
}
func (p plane) Pivot() int { return kdtree.Partition(p, kdtree.MedianOfRandoms(p, randoms)) }
func (p plane) Slice(start, end int) kdtree.SortSlicer {
	p.illuminaRecords = p.illuminaRecords[start:end]
	return p
}
func (p plane) Swap(i, j int) {
	p.illuminaRecords[i], p.illuminaRecords[j] = p.illuminaRecords[j], p.illuminaRecords[i]
}

// plane} OMIT

func main() {
	if len(os.Args) < 2 {
		fmt.Fprintln(os.Stderr, "missing input filename parameter")
		os.Exit(1)
	}

	meta := make(map[tileAddress]illuminaRecords)
	var discordant, total int
	{
		bf, err := boom.OpenBAM(os.Args[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "could not open file: %v\n", err)
			os.Exit(1)
		}
		names := bf.RefNames()

	load:
		for {
			var r [2]*boom.Record
			for i := range r {
				r[i], _, err = bf.Read()
				if err != nil {
					break load
				}
			}
			if r[0].Name() != r[1].Name() {
				fmt.Fprintln(os.Stderr, r[0].Name(), r[1].Name())
				fmt.Fprintln(os.Stderr, "name mismatch: file not sorted for name pairs")
				os.Exit(1)
			}
			total++

			const filterMask = boom.Unmapped | boom.MateUnmapped | boom.Secondary | boom.Duplicate | boom.ProperPair
			if r[0].Flags()&filterMask == 0 && r[1].Flags()&filterMask == 0 {
				discordant++
				m, err := newRecord(r, names)
				if err != nil {
					panic(err)
				}

				ta := tileAddress{
					FlowCell: m.FlowCell,
					Lane:     m.Lane,
					Tile:     m.Tile,
				}
				meta[ta] = append(meta[ta], m)
			}
		}
	}

	if len(meta) == 0 {
		fmt.Fprintln(os.Stderr, "no discordant read")
		os.Exit(0)
	}

	ts := make(map[tileAddress]*kdtree.Tree)
	// {build trees OMIT
	for ta, data := range meta {
		ts[ta] = kdtree.New(data, false)
	}
	// build trees} OMIT

	coincident := make([]int, len(offsets))
	{
		bf, err := boom.OpenBAM(os.Args[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "could not open file: %v\n", err)
			os.Exit(1)
		}
		names := bf.RefNames()

	search:
		for {
			var r [2]*boom.Record
			for i := range r {
				r[i], _, err = bf.Read()
				if err != nil {
					break search
				}
			}
			if r[0].Name() != r[1].Name() {
				fmt.Fprintln(os.Stderr, r[0].Name(), r[1].Name())
				panic("internal inconsistency: name mismatch")
			}

			// {concordant filter OMIT
			const (
				filterReq  = boom.ProperPair
				filterMask = boom.Unmapped | boom.MateUnmapped | boom.Secondary | boom.Duplicate | filterReq
			)
			// concordant filter} OMIT

			if r[0].Flags()&filterMask == filterReq && r[1].Flags()&filterMask == filterReq {
				q, err := newRecord(r, names)
				if err != nil {
					panic(err)
				}

				t, ok := ts[tileAddress{ // Get the relevant tree.
					FlowCell: q.FlowCell,
					Lane:     q.Lane,
					Tile:     q.Tile,
				}]
				if !ok { // We didn't have one, so there is no closest colony.
					continue
				}
				n, d := t.Nearest(q)
				if n == nil { // If there was a tree it must have a colony in it.
					panic("internal inconsistency: failed to find nearest")
				}
				nm := n.(*illuminaRecord)

				if nm.Metadata == q.Metadata { // We only stored discordant, only queried concordant.
					panic("internal inconsistency: discordant pair is concordant pair‽")
				}
				for i, off := range offsets {
					if overlap(q, nm, off.dist) {
						coincident[i]++
						fmt.Fprintf(os.Stderr, "@%d %0.fnm %+v -- %+v\n", off.dist, math.Sqrt(d), q, nm)
					}
				}
			}
		}
	}

	fmt.Printf("# %s\t%d\t%d\t%f\n",
		os.Args[1], total, discordant, float64(discordant)/float64(total),
	)
	for i, off := range offsets {
		fmt.Printf("%s\t%s\t%d\t%f\n",
			os.Args[1], off.label, coincident[i], float64(coincident[i])/float64(discordant),
		)
	}
}
