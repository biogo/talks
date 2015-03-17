// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math"
	"os"
	"sync"

	"github.com/biogo/boom"
	"github.com/biogo/illumina"
	"github.com/biogo/store/kdtree"
)

// Offset definition for overlap comparison.
type offset struct {
	dist  int
	label string
}

// The offsets we are interested in.
var offsets = []offset{
	{0, "Coincide"}, // The zeroth element is a special case - distributions are only kept for this offset.
	{1e2, "Adjacent"},
	{1e3, "At1k"},
	{1e4, "At10k"},
}

// This allows analysis of systems where there is some obliquity.
const (
	xunit = 37.5 // The width of a coordinate in nm.
	yunit = 37.5 // The height of a coordinate in nm.
)

// Index constants into the distance distribution table.
const (
	concord = iota
	discord
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
	A, B       mapping
	Concordant bool
	illumina.Metadata
}

// illuminaRecord} OMIT

// store is a string internment implementation.
type store map[string]string

// intern returns an interned version of the parameter.
func (is store) intern(s string) string {
	if s == "" {
		return ""
	}
	t, ok := is[s]
	if ok {
		return t
	}
	is[s] = s
	return s
}

var strings = make(store)

// newRecord returns an illumina record based on two boom.Records and a set of reference names.
func newRecord(r [2]*boom.Record, names []string) (*illuminaRecord, error) {
	m, err := illumina.Parse(boomIllumina{r[0]}) // They are a pair, so we only parse one.
	if err != nil {
		return nil, err
	}
	m.Instrument = strings.intern(m.Instrument)
	m.FlowCell = strings.intern(m.FlowCell)
	m.Multiplex.Tag = strings.intern(m.Multiplex.Tag)

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
		Concordant: r[0].Flags()&r[1].Flags()&boom.ProperPair != 0,
		Metadata:   m,
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

// buildTrees takes a map of illuminaRecords and returns a map of kdtrees, both keyed on
// tileAddress. Tree construction for each of the map elements is performed concurrently,
// allowing this part of the analysis to be performed in parallel.
func buildTrees(meta map[tileAddress]illuminaRecords) map[tileAddress]*kdtree.Tree {
	type tileTree struct {
		tile tileAddress
		tree *kdtree.Tree
	}

	r := make(chan tileTree, len(meta))

	var wg sync.WaitGroup
	wg.Add(len(meta))
	for ta, data := range meta {
		go func(ta tileAddress, data illuminaRecords) {
			defer wg.Done()
			r <- tileTree{
				tile: ta,
				tree: kdtree.New(data, false),
			}
		}(ta, data)
	}
	wg.Wait()
	close(r)

	ts := make(map[tileAddress]*kdtree.Tree)
	for t := range r {
		ts[t.tile] = t.tree
	}

	return ts
}

type dist []int

func newDist(r int) *dist {
	d := make(dist, r)
	return &d
}

func (d *dist) inc(i int) {
	s := *d
	switch {
	case i < len(s):
		s[i]++
		return
	case i == len(s):
		s = append(s, 1)
	case i < cap(s):
		s = s[:i+1]
		s[i] = 1
	case i >= cap(s):
		s = s[:cap(s)]
		s = append(s, make(dist, i+1-len(s))...)
		s[i] = 1
	}
	*d = s
}

var inf = math.Inf(1)

func main() {
	if len(os.Args) < 2 {
		fmt.Fprintln(os.Stderr, "missing input filename parameter")
		os.Exit(1)
	}

	meta := make(map[tileAddress]illuminaRecords)
	var discordant, concordant, mapped, total int
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

			const filterMask = boom.Unmapped | boom.MateUnmapped | boom.Secondary | boom.Duplicate
			if r[0].Flags()&filterMask == 0 && r[1].Flags()&filterMask == 0 {
				mapped++
				if r[0].Flags()&r[1].Flags()&boom.ProperPair != 0 {
					concordant++
				} else {
					discordant++
				}
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
		fmt.Fprintln(os.Stderr, "no mapped read")
		os.Exit(0)
	}

	// We do tree building concurrently with the larger case, c.f. discord-collision.go.
	ts := buildTrees(meta)

	var (
		concords = make([]int, len(offsets))
		discords = make([]int, len(offsets))

		allDist = make(map[tileAddress]*dist)
		dists   = make(map[tileAddress]*[2][2]*dist)
	)

	for ta := range ts {
		allDist[ta] = newDist(0)
		dists[ta] = &[2][2]*dist{
			{newDist(0), newDist(0)},
			{newDist(0), newDist(0)},
		}
	}

	{
		bf, err := boom.OpenBAM(os.Args[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "could not open file: %v\n", err)
			os.Exit(1)
		}
		names := bf.RefNames()

		nk := kdtree.NewNKeeper(2)
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

			const filterMask = boom.Unmapped | boom.MateUnmapped | boom.Secondary | boom.Duplicate
			if r[0].Flags()&filterMask == 0 && r[1].Flags()&filterMask == 0 {
				q, err := newRecord(r, names)
				if err != nil {
					panic(err)
				}

				ta := tileAddress{
					FlowCell: q.FlowCell,
					Lane:     q.Lane,
					Tile:     q.Tile,
				}
				t, ok := ts[ta]
				if !ok {
					continue
				}
				t.NearestSet(nk, q)
				if nk.Heap[0].Comparable == nil {
					panic("internal inconsistency: failed to find nearest")
				}
				if nk.Heap[1].Comparable == nil {
					// The second ComparableDist is the infinite distance marker,
					// so there was only one spot on the tile! We are it.
					continue
				}
				nm := nk.Heap[1].Comparable.(*illuminaRecord)
				d := int(math.Sqrt(nk.Heap[1].Dist))

				// Add the distance to the distribution for all mapped pairs, by tile.
				allDist[ta].inc(d)

				// Reset the keeper for the next query.
				nk.Heap = nk.Heap[:1]
				nk.Heap[0].Comparable = nil
				nk.Heap[0].Dist = inf

				if nm.Metadata == q.Metadata {
					panic("internal inconsistency: closest colony is second closest‽")
				}
				for i, off := range offsets {
					if overlap(q, nm, off.dist) {
						// Records concordant and discordant pairs as separate statistics.
						if i == 0 {
							// Keep distance distributions for the first offset class.
							if r[0].Flags()&r[1].Flags()&boom.ProperPair != 0 {
								concords[0]++
								if nm.Concordant {
									dists[ta][concord][concord].inc(d)
								} else {
									dists[ta][concord][discord].inc(d)
								}
							} else {
								discords[0]++
								if nm.Concordant {
									dists[ta][discord][concord].inc(d)
								} else {
									dists[ta][discord][discord].inc(d)
								}
							}
							fmt.Fprintf(os.Stderr, "%dnm %+v -- %+v\n", d, q, nm)
						} else {
							if r[0].Flags()&r[1].Flags()&boom.ProperPair != 0 {
								concords[i]++
							} else {
								discords[i]++
							}
						}
					}
				}
			}
		}
	}

	fmt.Printf("# %s\t%d\t%d\t%f\t%d\t%f\t%d\t%f\n",
		os.Args[1], total,
		mapped, float64(mapped)/float64(total),
		concordant, float64(concordant)/float64(total),
		discordant, float64(discordant)/float64(total),
	)

	for i, off := range offsets {
		fmt.Printf("%s\t%s\t%d\t%f\t%d\t%f\t%d\t%f\n",
			os.Args[1], off.label,
			concords[i]+discords[i], float64(concords[i]+discords[i])/float64(mapped),
			concords[i], float64(concords[i])/float64(concordant),
			discords[i], float64(discords[i])/float64(discordant),
		)
	}

	for ta := range ts {
		for _, ds := range []struct {
			label string
			dist
		}{
			{"All", *allDist[ta]},
			{"Concord Concord", *dists[ta][concord][concord]},
			{"Concord Discord", *dists[ta][concord][discord]},
			{"Discord Concord", *dists[ta][discord][concord]},
			{"Discord Discord", *dists[ta][discord][discord]},
		} {
			for d, c := range ds.dist {
				if c != 0 {
					fmt.Printf("%s\t%s.%d.%d\t%s\t%d\t%d\n",
						os.Args[1], ta.FlowCell, ta.Lane, ta.Tile, ds.label, d, c,
					)
				}
			}
		}
	}
}
