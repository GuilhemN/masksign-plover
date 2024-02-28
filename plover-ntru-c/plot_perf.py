import re

pgf_data = ""
pgf_data += "d  kg_ms  kg_kB  sg_ms  sg_kB  vf_ms  vf_kB\n"

tex = ""
tex += "\\textbf{Variant} & & \\textbf{Keygen} & & & \\textbf{Sign} & & & \\textbf{Verify} & \\\\\n"
tex += "$\kappa-d$ & ms & Mclk & stack & ms & Mclk & stack & ms & Mclk & stack\\\\\n"
tex += "\\hline\n"

# first parse data
for d in [1, 2, 4, 8, 16, 32]:
    f = open(f"bench_PLOVER_128_{d}.txt", "r")
    content = ''.join(f.readlines())
    # print(content)

    # extract timings
    m = re.search(
f"""Plover-128-{d}\\s+KeyGen\\(\\)\\s+[0-9]+:\\s+(?P<keygen_ms>[0-9.]+) ms\\s+(?P<keygen_mcyc>[0-9.]+) Mcyc
Plover-128-{d}\\s+Sign\(\)\\s+[0-9]+:\\s+(?P<sign_ms>[0-9.]+) ms\\s+(?P<sign_mcyc>[0-9.]+) Mcyc
Plover-128-{d}\\s+Verify\(\)\\s+[0-9]+:\\s+(?P<verify_ms>[0-9.]+) ms\\s+(?P<verify_mcyc>[0-9.]+) Mcyc""", content)
    keygen_ms, keygen_mcyc = m.group("keygen_ms"), m.group("keygen_mcyc")
    sign_ms, sign_mcyc = m.group("sign_ms"), m.group("sign_mcyc")
    verify_ms, verify_mcyc = m.group("verify_ms"), m.group("verify_mcyc")

    # extract stack size
    m = re.search(
f"""plover_core.c:[0-9:A-Za-z_]+_keygen\\s+(?P<keygen_stack>[0-9]+)\\s+static
plover_core.c:[0-9:A-Za-z_]+_sign\\s+(?P<sign_stack>[0-9]+)\\s+static
plover_core.c:[0-9:A-Za-z_]+_verify\\s+(?P<verify_stack>[0-9]+)\\s+static""", content)
    keygen_stack = m.group("keygen_stack")
    sign_stack = m.group("sign_stack")
    verify_stack = m.group("verify_stack")

    f.close()

    pgf_data += f"{d:2} {keygen_ms:6} {int(keygen_stack)//1000:5}  {sign_ms:6}  {int(sign_stack)//1000:4}  {verify_ms:5}     {int(verify_stack)//1000}\n"

    if d != 1:
        verify_ms = verify_mcyc = verify_stack = "="
    tex += f"128-{d} & {keygen_ms} & {keygen_mcyc} & {keygen_stack} & {sign_ms} & {sign_mcyc} & {sign_stack} & {verify_ms} & {verify_mcyc} & {verify_stack} \\\\\n"

# print(tex)
print(pgf_data)