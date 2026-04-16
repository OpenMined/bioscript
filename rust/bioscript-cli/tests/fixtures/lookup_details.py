VARIANT = bioscript.variant(
    rsid="rs73885319",
    grch38="22:36265860-36265860",
    ref="A",
    alt="G",
    kind="snp",
)


def main():
    genotypes = bioscript.load_genotypes(input_file)
    details = genotypes.lookup_variant_details(VARIANT)
    print(details)


if __name__ == "__main__":
    main()
