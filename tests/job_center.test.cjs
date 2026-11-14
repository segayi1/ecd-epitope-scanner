const assert = require("node:assert/strict");
const test = require("node:test");

const {
  jobLabel,
  partitionJobs,
  progressPercent,
} = require("../web/job_center.js");

test("active jobs are ordered before recent terminal jobs", () => {
  const jobs = [
    { id: "complete", status: "COMPLETE" },
    { id: "queued", status: "QUEUED" },
    { id: "running", status: "RUNNING" },
    { id: "waiting", status: "WAITING_FOR_AF3_UPLOAD" },
    { id: "failed", status: "FAILED" },
  ];

  const result = partitionJobs(jobs, 1);

  assert.deepEqual(result.active.map((job) => job.id), [
    "running",
    "queued",
    "waiting",
  ]);
  assert.deepEqual(result.recent.map((job) => job.id), ["complete"]);
});

test("task labels and progress tolerate missing or odd values", () => {
  assert.equal(jobLabel({ id: "job-1", identifier: "P16410" }), "P16410");
  assert.equal(jobLabel({ id: "job-1" }), "job-1");
  assert.equal(progressPercent(0.626), 63);
  assert.equal(progressPercent("not-a-number"), 0);
  assert.equal(progressPercent(1.5), 100);
});
