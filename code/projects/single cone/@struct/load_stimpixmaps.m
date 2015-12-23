function maps = load_stimpixmaps(datarun)

parsed = parse_rrs_prefix(datarun);
maps = load_stimpixmaps(parsed.piece_fullname);