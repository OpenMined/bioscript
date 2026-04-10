base = "abc"
result = base.replace("a", "z" * 1000000)
print(len(result))
