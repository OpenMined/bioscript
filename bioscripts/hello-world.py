def add(a, b):
    return a + b


def main():
    total = add(2, 3)
    print("hello from bioscript")
    print("2 + 3 =", total)

    source = read_text("bioscripts/input.txt").strip()
    print("loaded:", source)

    write_text(
        "bioscripts/output/hello-world.txt",
        "hello from bioscript\n"
        + "2 + 3 = "
        + str(total)
        + "\n"
        + "loaded: "
        + source
        + "\n",
    )


if __name__ == "__main__":
    main()
