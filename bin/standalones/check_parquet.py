import duckdb

parquet_file = "/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/IPS/asgard_ips_results_feb5.parquet"

con = duckdb.connect()

print(
    con.execute(
       f"SELECT * FROM read_parquet('{parquet_file}') LIMIT 5"
    ).fetchdf(),

    con.execute(
        f"DESCRIBE SELECT * FROM read_parquet('{parquet_file}')"
    ).fetchdf()
)