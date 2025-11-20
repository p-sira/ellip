/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

let camera = document.getElementById("plotly-html-element")._fullLayout.scene.camera;
let center = camera.center;
let eye = camera.eye;
console.log([center.x, center.y, center.z]);
console.log([eye.x, eye.y, eye.z]);
