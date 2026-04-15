document.addEventListener("DOMContentLoaded", () => {
  const form = document.getElementById("motiffinder-form");
  const status = document.getElementById("status");
  const output = document.getElementById("output");

  form.addEventListener("submit", async (e) => {
    e.preventDefault();
    status.textContent = "Running MotifFinder…";
    output.innerHTML = "";

    const genbank = document.getElementById("genbank").files[0];
    const motifs = document.getElementById("motifs").files[0];

    if (!genbank || !motifs) {
      status.textContent = "Please select both files.";
      return;
    }

    const data = new FormData();
    data.append("genbank", genbank);
    data.append("motifs", motifs);

    try {
      const res = await fetch("/api/run_motiffinder", {
        method: "POST",
        body: data,
      });

      const json = await res.json();
      status.textContent = `Found ${json.count} motif hits`;

      json.results.forEach((r) => {
        const line = document.createElement("div");
        line.textContent =
          `${r.strand} ${r.start}-${r.end} ` +
          `motif=${r.motif} features=${r.features.join(",")}`;
        output.appendChild(line);
      });
    } catch (err) {
      status.textContent = "Error running MotifFinder.";
      console.error(err);
    }
  });
});