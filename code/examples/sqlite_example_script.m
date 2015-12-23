dbid = mksqlite('open', snle_sqlite_db_path) % If only one DB is open, DBID should be 1
tables = mksqlite(dbid, 'show tables'); % Includes tables and views
temporal_pieces = mksqlite(dbid, 'SELECT piece FROM eccentricity WHERE piece_ecc > 5 AND clock < 6'); % If DBID == 1, this argument can be left out