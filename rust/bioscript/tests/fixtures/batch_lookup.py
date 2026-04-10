VARIANT_A = bioscript.variant(
    rsid="rs73885319",
    grch38="22:36265860-36265860",
    ref="A",
    alt="G",
    kind="snp",
)

VARIANT_B = bioscript.variant(
    rsid="rs60910145",
    grch38="22:36265988-36265988",
    ref="T",
    alt="G",
    kind="snp",
)

VARIANT_C = bioscript.variant(
    rsid=["rs71785313", "rs1317778148", "rs143830837"],
    grch38="22:36266000-36266005",
    ref="I",
    alt="D",
    kind="deletion",
    deletion_length=6,
)

QUERY_PLAN = bioscript.query_plan([
    VARIANT_A,
    VARIANT_B,
    VARIANT_C,
])


def main():
    genotypes = bioscript.load_genotypes(input_file)
    site_a, site_b, site_c = genotypes.lookup_variants(QUERY_PLAN)
    print(site_a)
    print(site_b)
    print(site_c)


if __name__ == "__main__":
    main()
