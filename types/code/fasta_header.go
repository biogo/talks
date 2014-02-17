func (r *Reader) header(line []byte) (seqio.SequenceAppender, error) {
	s := r.t.Clone().(seqio.SequenceAppender)
	fieldMark := bytes.IndexAny(line, " \t")
	var err error
	if fieldMark < 0 {
		err = s.SetName(string(line[len(r.IDPrefix):]))
		return s, err
	} else {
		err = s.SetName(string(line[len(r.IDPrefix):fieldMark]))
		_err := s.SetDescription(string(line[fieldMark+1:])) // HL
		if err != nil || _err != nil {
			switch {
			case err == _err:
				return s, err
			case err != nil && _err != nil:
				return s, fmt.Errorf("fasta: multiple errors: name: %s, desc:%s", err, _err)
			case err != nil:
				return s, err
			case _err != nil:
				return s, _err
			}
		}
	}

	return s, nil
}
