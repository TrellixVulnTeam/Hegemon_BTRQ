import os

def conf_file(hegemon_base: str, output: str) -> None:
    export = os.path.join(output, hegemon_base + "-explore.txt")
    with open(export, "w") as file_out:
        file_out.write("[]\n")
        file_out.write("name=\n")

        names = ["expr", "index", "survival", "indexHeader", "info"]
        types = ["expr", "idx", "survival", "ih", "info"]
        for name, type in zip(names, types):
            file = f"{hegemon_base}-{type}.txt"
            filepath = os.path.join(output, file)
            file_out.write(f"{name}={filepath}\n")

        file_out.write("key=")
    print(f"Conf file made at {output}")