// Copyright ©2013 The bíogo.talks Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//+build ignore

package main

import (
	"fmt"

	"github.com/biogo/illumina"

	"github.com/davecgh/go-spew/spew"
)

type Read struct {
	name string
	desc string
}

func (t Read) Name() string        { return t.name }
func (t Read) Description() string { return t.desc }

func main() {
	reads := []Read{
		{"HWUSI-EAS100R:6:73:941:1973#ATCACG/1", ""},
		{"EAS139:136:FC706VJ:2:2104:15343:197393", "1:Y:18:ATCACG"},
	}

	for _, r := range reads {
		m, err := illumina.Parse(r) // HL
		if err != nil {
			fmt.Println(err)
		} else {
			spew.Dump(m)
		}
	}
}
