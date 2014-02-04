// Copyright ©2013 The bíogo.talks Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"code.google.com/p/biogo.store/kdtree"

	"fmt"
	"math"
)

var wpData = kdtree.Points{{2, 3}, {5, 4}, {9, 6}, {4, 7}, {8, 1}, {7, 2}}

func main() {
	t := kdtree.New(wpData, false)
	q := kdtree.Point{8, 7}
	p, d := t.Nearest(q)
	fmt.Printf("%v is closest point to %v, d=%f\n", p, q, math.Sqrt(d))
}
