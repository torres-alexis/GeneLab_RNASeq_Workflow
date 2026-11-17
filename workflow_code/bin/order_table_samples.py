#!/usr/bin/env python3

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description="Order table sample columns by the runsheet")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--runsheet")
    args = parser.parse_args()

    with open(args.samples) as handle:
        samples = handle.read().splitlines()
    aliases = {}
    if args.runsheet:
        with open(args.runsheet, newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("Has Tech Reps", "").strip().upper() == "TRUE":
                    aliases[row["Sample Name"]] = row.get("Source Name", "")

    with open(args.input, newline="") as source, open(args.output, "w", newline="") as output:
        reader = csv.reader(source)
        writer = csv.writer(output)
        header = next(reader)
        uncollapsed = {aliases[name] for name in header[1:] if name in aliases}
        ordered_names = []
        for name in samples:
            ordered_names.append(name)
            source_name = aliases.get(name)
            if source_name and source_name not in uncollapsed:
                ordered_names.append(source_name)
        sample_order = {name: i for i, name in enumerate(dict.fromkeys(ordered_names))}
        slots = [i for i, name in enumerate(header) if i > 0 and name in sample_order]
        ordered = sorted(slots, key=lambda i: sample_order[header[i]])
        columns = list(range(len(header)))
        for slot, index in zip(slots, ordered):
            columns[slot] = index
        writer.writerow([header[i] for i in columns])
        for row in reader:
            if not row:
                writer.writerow(row)
                continue
            row += [""] * (len(header) - len(row))
            writer.writerow([row[i] for i in columns] + row[len(header):])


if __name__ == "__main__":
    main()
