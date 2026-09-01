/**
 * Shared file upload and drop-zone handlers for all SyToGen form pages.
 * Used by codonbias.html, motiffinder.html, and sytogen.html.
 */

/**
 * Attach click and drag-drop handlers to a file input via its drop zone.
 * @param {string} dzId - ID of the drop zone container
 * @param {string} inputId - ID of the hidden file input
 * @param {string} dpId - ID of the drop-primary text element (optional, for label updates)
 */
function attachDropZone(dzId, inputId, dpId) {
  const dz = document.getElementById(dzId);
  const input = document.getElementById(inputId);
  const dp = dpId ? document.getElementById(dpId) : null;

  if (!dz || !input) return;

  // Click to open file picker
  dz.addEventListener("click", () => input.click());

  // File selected via picker
  input.addEventListener("change", () => {
    if (input.files.length && dp) {
      dp.textContent = input.files[0].name;
    }
  });

  // Drag over
  dz.addEventListener("dragover", (event) => {
    event.preventDefault();
    event.stopPropagation();
    dz.classList.add("drag-active");
  });

  // Drag leave
  dz.addEventListener("dragleave", (event) => {
    event.preventDefault();
    event.stopPropagation();
    dz.classList.remove("drag-active");
  });

  // Drop
  dz.addEventListener("drop", (event) => {
    event.preventDefault();
    event.stopPropagation();
    dz.classList.remove("drag-active");

    const files = event.dataTransfer.files;
    if (!files || !files.length) return;

    input.files = files;
    input.dispatchEvent(new Event("change", { bubbles: true }));
  });
}

/**
 * Update a label element when a file is selected.
 * @param {string} inputId - ID of the file input
 * @param {string} labelId - ID of the label element to update
 */
function bindLabel(inputId, labelId) {
  const input = document.getElementById(inputId);
  const label = document.getElementById(labelId);

  if (!input || !label) return;

  input.addEventListener("change", () => {
    if (input.files.length) {
      label.textContent = input.files[0].name;
    }
  });
}

/**
 * Attach drag-and-drop handlers to a zone element for file selection.
 * @param {string} zoneId - ID of the drop zone container
 * @param {string} inputId - ID of the hidden file input
 */
function bindDropZone(zoneId, inputId) {
  const zone = document.getElementById(zoneId);
  const input = document.getElementById(inputId);

  if (!zone || !input) return;

  ["dragenter", "dragover"].forEach((evt) => {
    zone.addEventListener(evt, (event) => {
      event.preventDefault();
      event.stopPropagation();
      zone.classList.add("drag-active");
    });
  });

  ["dragleave", "dragend"].forEach((evt) => {
    zone.addEventListener(evt, (event) => {
      event.preventDefault();
      event.stopPropagation();
      zone.classList.remove("drag-active");
    });
  });

  zone.addEventListener("drop", (event) => {
    event.preventDefault();
    event.stopPropagation();
    zone.classList.remove("drag-active");

    const files = event.dataTransfer.files;
    if (!files || !files.length) return;

    input.files = files;
    input.dispatchEvent(new Event("change", { bubbles: true }));
  });
}
