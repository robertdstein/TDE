from tabulate import tabulate
import os

def run(table, headers, title, savename):

    print savename

    path = "tex_table/" + savename + ".tex"
    with open(path, "w") as f:

        # variables = ["Name", "Redshift", "RA", "Dec", "Max Date"]

        f.write("\documentclass[]{article} \n")
        f.write("\usepackage{siunitx} \n")
        f.write(r"\begin{document}" + "\n")
        f.write(r"\begin{tabular}{ |p{3.5cm}||" + (len(headers) - 1) * \
                "p{4.5cm}|" + "} \n")
        f.write("\hline \n")
        f.write("\multicolumn{" + str(len(headers)) + "}{|c|}{" + title +
                r"} \\ " + "\n")

        f.write("\hline \n")
        f.write("&".join([x for x in headers]) + r" \\" + " \n")
        f.write("\hline \n")
        for row in table:
            line = ""

            for entry in row:
                if isinstance(entry, float):
                    line += r" \num[round-precision=2, round-mode=figures, " + \
                            r"scientific-notation=true]{" + str(entry) + "} "
                else:
                    line += entry
                line += " &"

            line = line[:-1] + r"\\" + " \n"

            f.write(line)

        f.write("\hline \n")
        f.write("\end{tabular} \n")
        f.write("\end{document} \n")


    # os.system("pdflatex " + path)


