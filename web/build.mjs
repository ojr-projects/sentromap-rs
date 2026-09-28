// Bundle src/main.ts into dist/app.js and copy static files. `--watch` rebuilds on change.
import * as esbuild from "esbuild";
import { cpSync, mkdirSync } from "node:fs";

mkdirSync("dist", { recursive: true });
const copy = () => {
  cpSync("index.html", "dist/index.html");
  cpSync("src/style.css", "dist/style.css");
};
const opts = {
  entryPoints: ["src/main.ts"],
  bundle: true,
  format: "esm",
  target: "es2022",
  outfile: "dist/app.js",
  sourcemap: true,
  minify: !process.argv.includes("--watch"),
  plugins: [{ name: "copy", setup: (b) => b.onEnd(copy) }],
};
if (process.argv.includes("--watch")) {
  const ctx = await esbuild.context(opts);
  await ctx.watch();
  console.log("watching…");
} else {
  await esbuild.build(opts);
}
