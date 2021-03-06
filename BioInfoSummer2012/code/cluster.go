// Copyright ©2012 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"strings"

	"github.com/biogo/cluster/kmeans" // @090d6
)

const scaling = 30

type Feature struct {
	id    string
	start int
	end   int
}

func (f *Feature) Start() int { return f.start }
func (f *Feature) End() int   { return f.end }
func (f *Feature) Len() int   { return f.end - f.start }
func (f *Feature) String() string {
	return fmt.Sprintf("%2s %s%s",
		f.id,
		strings.Repeat(" ", f.Start()/scaling),
		strings.Repeat("-", f.Len()/scaling),
	)
}

type Features []*Feature

func (f Features) Len() int { return len(f) }
func (f Features) Values(i int) []float64 {
	return []float64{float64(f[i].Start()), float64(f[i].End())}
}

var feats = []*Feature{
	{id: "0", start: 1, end: 1700},
	{id: "1", start: 3, end: 610},
	{id: "2", start: 2, end: 750},
	{id: "3", start: 650, end: 900},
	{id: "4", start: 2, end: 1700},
	{id: "5", start: 700, end: 950},
	{id: "6", start: 2, end: 605},
	{id: "7", start: 1000, end: 1700},
	{id: "8", start: 950, end: 1712},
	{id: "9", start: 1, end: 600},
	{id: "10", start: 1000, end: 1650},
}

// Cluster feat.Features on the basis of location where:
//  epsilon is allowable error, and
//  effort is number of attempts to achieve error < epsilon for any k.
func ClusterFeatures(f []*Feature, epsilon float64, effort int) *kmeans.Kmeans {
	km, err := kmeans.New(Features(f))
	if err != nil {
		panic(err)
	}

	values := km.Values()
	cut := make([]float64, len(values))
	for i, v := range values {
		v := v.V()
		l := epsilon * (v[1] - v[0])
		cut[i] = l * l
	}

	for k := 1; k <= len(f); k++ {
	ATTEMPT:
		for attempt := 0; attempt < effort; attempt++ {
			km.Seed(k)
			km.Cluster()
			centers := km.Centers()
			for i, v := range values {
				cv := centers[v.Cluster()].V()
				vv := v.V()
				dx, dy := cv[0]-vv[0], cv[1]-vv[1]
				ok := dx*dx+dy*dy < cut[i]
				if !ok {
					continue ATTEMPT
				}
			}
			return km
		}
	}

	panic("cannot reach")
}

func main() {
	km := ClusterFeatures(feats, 0.15, 5)
	for ci, c := range km.Centers() {
		fmt.Printf("Cluster %d:\n", ci)
		for _, i := range c.Cluster() {
			fmt.Println(feats[i])
		}
		fmt.Println()
	}

	var within float64
	for _, ss := range km.Within() {
		within += ss
	}
	fmt.Printf("betweenSS / totalSS = %.6f\n", 1-(within/km.Total()))
}
