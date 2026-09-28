// Headless screenshot for UI smoke checks: node shot.mjs URL OUT.png [width height] [waitMs]
import { chromium } from "playwright";

const [url, out, w = "1600", h = "1000", wait = "2500"] = process.argv.slice(2);
const browser = await chromium.launch();
const page = await browser.newPage({ viewport: { width: Number(w), height: Number(h) } });
const errors = [];
page.on("console", (m) => m.type() === "error" && errors.push(m.text()));
page.on("pageerror", (e) => errors.push(String(e)));
await page.goto(url);
await page.waitForTimeout(Number(wait));
// Optional pointer actions: HOVER="x,y" moves the mouse there; CLICK="x,y" clicks.
for (const [env, act] of [["HOVER", "move"], ["CLICK", "click"]]) {
  if (!process.env[env]) continue;
  const [x, y] = process.env[env].split(",").map(Number);
  if (act === "move") await page.mouse.move(x, y);
  else await page.mouse.click(x, y);
  await page.waitForTimeout(800);
}
await page.screenshot({ path: out });
if (errors.length) console.log("page errors:\n" + errors.join("\n"));
await browser.close();
