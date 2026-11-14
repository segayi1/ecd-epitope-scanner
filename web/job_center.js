(() => {
  "use strict";

  const ACTIVE_RANK = new Map([
    ["RUNNING", 0],
    ["CANCELLING", 0],
    ["QUEUED", 1],
    ["WAITING_FOR_AF3_UPLOAD", 2],
  ]);
  const RECENT_LIMIT = 8;

  function partitionJobs(jobs, recentLimit = RECENT_LIMIT) {
    const source = Array.isArray(jobs) ? jobs : [];
    const active = source
      .filter((job) => ACTIVE_RANK.has(job.status))
      .map((job, index) => ({ job, index }))
      .sort((left, right) =>
        ACTIVE_RANK.get(left.job.status) - ACTIVE_RANK.get(right.job.status) ||
        left.index - right.index
      )
      .map(({ job }) => job);
    const recent = source
      .filter((job) => !ACTIVE_RANK.has(job.status))
      .slice(0, recentLimit);
    return { active, recent };
  }

  function progressPercent(value) {
    const number = Number(value);
    if (!Number.isFinite(number)) return 0;
    return Math.max(0, Math.min(100, Math.round(number * 100)));
  }

  function jobLabel(job) {
    return job.accession || job.identifier || job.id;
  }

  if (typeof document === "undefined") {
    if (typeof module !== "undefined") {
      module.exports = { partitionJobs, progressPercent, jobLabel };
    }
    return;
  }

  const COPY = {
    zh: {
      tasks: "任务",
      title: "分析任务",
      intro: "运行中、排队中和最近结束的分析会自动刷新。",
      active: "运行中与排队",
      recent: "最近结束",
      emptyActive: "当前没有运行或排队中的分析。",
      emptyRecent: "暂无最近完成的分析。",
      open: "打开",
      progress: "进度",
      updated: "更新",
      loadError: "任务列表暂时无法读取。",
      cancel: "终止任务",
      cancellingMessage: "已请求终止，等待当前计算到达可停止的位置。",
      cancelError: "终止请求失败，请刷新后重试。",
      close: "关闭",
      status: {
        RUNNING: "运行中",
        CANCELLING: "正在终止",
        CANCELLED: "已终止",
        QUEUED: "排队中",
        WAITING_FOR_AF3_UPLOAD: "等待 AF3",
        COMPLETE: "完成",
        FAILED: "失败",
      },
    },
    en: {
      tasks: "Tasks",
      title: "Analysis tasks",
      intro: "Running, queued and recently finished analyses refresh automatically.",
      active: "Running and queued",
      recent: "Recently finished",
      emptyActive: "No analysis is running or queued.",
      emptyRecent: "No recent finished analysis.",
      open: "Open",
      progress: "Progress",
      updated: "Updated",
      loadError: "The task list is temporarily unavailable.",
      cancel: "Stop task",
      cancellingMessage: "Stop requested; waiting for the next cancellation checkpoint.",
      cancelError: "Could not stop the task; refresh and try again.",
      close: "Close",
      status: {
        RUNNING: "Running",
        CANCELLING: "Stopping",
        CANCELLED: "Stopped",
        QUEUED: "Queued",
        WAITING_FOR_AF3_UPLOAD: "Waiting for AF3",
        COMPLETE: "Complete",
        FAILED: "Failed",
      },
    },
  };

  const STYLE = `
    .job-center-trigger {
      position: relative;
      display: inline-flex;
      align-items: center;
      gap: 7px;
      min-height: 34px;
      padding: 7px 10px;
      border: 1px solid var(--line-strong, #c8d5d9);
      border-radius: 8px;
      background: #fff;
      color: var(--ink, #172b35);
      font: inherit;
      font-weight: 650;
      cursor: pointer;
      white-space: nowrap;
    }
    .job-center-trigger:hover,
    .job-center-trigger:focus-visible {
      border-color: var(--teal, #127c80);
      outline: none;
      box-shadow: 0 0 0 3px #127c801a;
    }
    .job-center-trigger.has-active {
      border-color: var(--teal, #127c80);
      background: #edf8f6;
      color: var(--teal-dark, #0b5c60);
    }
    .job-center-count {
      display: grid;
      min-width: 20px;
      height: 20px;
      place-items: center;
      padding: 0 5px;
      border-radius: 999px;
      background: var(--teal, #127c80);
      color: #fff;
      font-size: 11px;
      line-height: 1;
    }
    .job-center-count[hidden] { display: none; }
    .job-center-dialog {
      width: min(680px, calc(100vw - 28px));
      max-height: min(760px, calc(100vh - 32px));
      padding: 0;
      overflow: hidden;
      border: 1px solid #cbd7db;
      border-radius: 16px;
      background: #f7faf9;
      color: var(--ink, #172b35);
      box-shadow: 0 24px 70px #10262e42;
    }
    .job-center-dialog::backdrop {
      background: #10262e66;
      backdrop-filter: blur(2px);
    }
    .job-center-head {
      display: flex;
      align-items: flex-start;
      justify-content: space-between;
      gap: 20px;
      padding: 20px 22px 15px;
      border-bottom: 1px solid #dce5e8;
      background: #fff;
    }
    .job-center-head h2 {
      margin: 0;
      font-size: 1.25rem;
      letter-spacing: -.018em;
    }
    .job-center-head p {
      margin: 5px 0 0;
      color: #657780;
      font-size: 12px;
    }
    .job-center-close {
      width: 34px;
      height: 34px;
      flex: 0 0 auto;
      border: 1px solid #d5e0e3;
      border-radius: 9px;
      background: #fff;
      color: #52686f;
      font: 20px/1 sans-serif;
      cursor: pointer;
    }
    .job-center-body {
      max-height: calc(100vh - 170px);
      overflow: auto;
      padding: 18px 22px 22px;
    }
    .job-center-section + .job-center-section { margin-top: 22px; }
    .job-center-section h3 {
      margin: 0 0 9px;
      color: #52686f;
      font-size: 12px;
      letter-spacing: .08em;
      text-transform: uppercase;
    }
    .job-center-list {
      display: grid;
      gap: 8px;
      margin: 0;
      padding: 0;
      list-style: none;
    }
    .job-center-item {
      display: grid;
      grid-template-columns: minmax(0, 1fr) auto;
      gap: 9px 16px;
      width: 100%;
      padding: 13px 14px;
      border: 1px solid #d8e2e5;
      border-radius: 11px;
      background: #fff;
      color: inherit;
      font: inherit;
      text-align: left;
      cursor: pointer;
    }
    .job-center-item:hover,
    .job-center-item:focus-visible {
      border-color: #127c80;
      outline: none;
      box-shadow: 0 0 0 3px #127c8017;
    }
    .job-center-item-main {
      display: grid;
      min-width: 0;
      gap: 3px;
    }
    .job-center-item-title {
      overflow: hidden;
      font-weight: 700;
      text-overflow: ellipsis;
      white-space: nowrap;
    }
    .job-center-meta,
    .job-center-message {
      overflow: hidden;
      color: #657780;
      font-size: 12px;
      text-overflow: ellipsis;
      white-space: nowrap;
    }
    .job-center-side {
      display: grid;
      align-content: start;
      justify-items: end;
      gap: 4px;
    }
    .job-center-status {
      padding: 3px 7px;
      border-radius: 999px;
      background: #edf1f2;
      color: #52616a;
      font-size: 11px;
      font-weight: 700;
    }
    .job-center-status.running,
    .job-center-status.queued,
    .job-center-status.waiting {
      background: #dff2ee;
      color: #0b666a;
    }
    .job-center-status.failed {
      background: #fde9e6;
      color: #963a31;
    }
    .job-center-track {
      grid-column: 1 / -1;
      height: 5px;
      overflow: hidden;
      border-radius: 999px;
      background: #e7edef;
    }
    .job-center-bar {
      display: block;
      height: 100%;
      border-radius: inherit;
      background: #127c80;
      transition: width .25s ease;
    }
    .job-center-empty,
    .job-center-error {
      margin: 0;
      padding: 14px;
      border: 1px dashed #cbd7db;
      border-radius: 10px;
      color: #657780;
      font-size: 13px;
    }
    .job-center-error { color: #963a31; }
    .job-center-cancel { margin: 6px 0 2px; padding: 6px 10px; border: 1px solid #b65a50; border-radius: 7px; background: #fff; color: #963a31; cursor: pointer; font: inherit; }
    @media (max-width: 620px) {
      .job-center-dialog { width: calc(100vw - 16px); max-height: calc(100vh - 16px); }
      .job-center-head { padding: 16px; }
      .job-center-body { max-height: calc(100vh - 135px); padding: 14px 16px 16px; }
      .job-center-item { gap: 7px 10px; padding: 12px; }
    }
  `;

  const state = { jobs: [], error: "", loading: false };

  function language() {
    try {
      return localStorage.getItem("ecd.lang") === "en" ? "en" : "zh";
    } catch {
      return "zh";
    }
  }

  function copy() {
    return COPY[language()];
  }

  function statusText(status) {
    return copy().status[status] || status;
  }

  function statusClass(status) {
    return {
      RUNNING: "running",
      QUEUED: "queued",
      WAITING_FOR_AF3_UPLOAD: "waiting",
      FAILED: "failed",
    }[status] || "";
  }

  function updatedText(value) {
    if (!value) return "";
    const date = new Date(value);
    if (Number.isNaN(date.getTime())) return "";
    return new Intl.DateTimeFormat(language() === "en" ? "en" : "zh-CN", {
      month: "short",
      day: "numeric",
      hour: "2-digit",
      minute: "2-digit",
    }).format(date);
  }

  function build() {
    if (document.getElementById("global-job-center")) return;

    const style = document.createElement("style");
    style.textContent = STYLE;
    document.head.append(style);

    const trigger = document.createElement("button");
    trigger.id = "global-job-center-trigger";
    trigger.className = "job-center-trigger";
    trigger.type = "button";
    trigger.setAttribute("aria-haspopup", "dialog");
    trigger.setAttribute("aria-controls", "global-job-center");
    trigger.innerHTML = '<span data-job-center-label></span><span class="job-center-count" hidden></span>';

    const host = document.querySelector(".app-actions") || document.querySelector(".site-links");
    if (host) {
      const runButton = host.querySelector("#open-run");
      host.insertBefore(trigger, runButton || null);
    }

    const dialog = document.createElement("dialog");
    dialog.id = "global-job-center";
    dialog.className = "job-center-dialog";
    dialog.innerHTML = `
      <div class="job-center-head">
        <div><h2 data-job-center-title></h2><p data-job-center-intro></p></div>
        <button class="job-center-close" type="button" data-job-center-close aria-label="">×</button>
      </div>
      <div class="job-center-body">
        <section class="job-center-section">
          <h3 data-job-center-active></h3>
          <div data-job-center-active-list></div>
        </section>
        <section class="job-center-section">
          <h3 data-job-center-recent></h3>
          <div data-job-center-recent-list></div>
        </section>
        <p class="job-center-error" data-job-center-error hidden></p>
      </div>`;
    document.body.append(dialog);

    trigger.addEventListener("click", () => {
      render();
      if (typeof dialog.showModal === "function") dialog.showModal();
      else dialog.setAttribute("open", "");
      loadJobs();
    });
    dialog.querySelector("[data-job-center-close]").addEventListener("click", () => dialog.close());
    dialog.addEventListener("click", (event) => {
      if (event.target === dialog) dialog.close();
    });
    window.addEventListener("storage", (event) => {
      if (event.key === "ecd.lang") render();
    });
    document.addEventListener("visibilitychange", () => {
      if (document.visibilityState === "visible") loadJobs();
    });
    const languageButton = document.getElementById("btn-lang");
    if (languageButton) languageButton.addEventListener("click", () => setTimeout(render, 0));

    window.setInterval(loadJobs, 2500);
    loadJobs();
  }

  function renderList(container, jobs, emptyText) {
    container.replaceChildren();
    if (!jobs.length) {
      const empty = document.createElement("p");
      empty.className = "job-center-empty";
      empty.textContent = emptyText;
      container.append(empty);
      return;
    }
    const list = document.createElement("ul");
    list.className = "job-center-list";
    jobs.forEach((job) => {
      const item = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.className = "job-center-item";
      button.title = `${job.id} · ${jobLabel(job)}`;
      button.addEventListener("click", () => {
        location.assign(`/?job=${encodeURIComponent(job.id)}`);
      });

      const main = document.createElement("span");
      main.className = "job-center-item-main";
      const title = document.createElement("span");
      title.className = "job-center-item-title";
      title.textContent = jobLabel(job);
      const meta = document.createElement("span");
      meta.className = "job-center-meta";
      meta.textContent = `${job.id} · ${copy().updated} ${updatedText(job.updated_at) || "-"}`;
      const message = document.createElement("span");
      message.className = "job-center-message";
      message.textContent = job.status === "CANCELLING" ? copy().cancellingMessage : job.error || job.message || job.stage || "";
      main.append(title, meta, message);

      const side = document.createElement("span");
      side.className = "job-center-side";
      const status = document.createElement("span");
      status.className = `job-center-status ${statusClass(job.status)}`.trim();
      status.textContent = statusText(job.status);
      const open = document.createElement("span");
      open.className = "job-center-meta";
      open.textContent = `${copy().open} →`;
      side.append(status, open);

      button.append(main, side);
      if (ACTIVE_RANK.has(job.status)) {
        const percent = progressPercent(job.progress);
        const track = document.createElement("span");
        track.className = "job-center-track";
        track.setAttribute("role", "progressbar");
        track.setAttribute("aria-label", `${copy().progress} ${percent}%`);
        track.setAttribute("aria-valuemin", "0");
        track.setAttribute("aria-valuemax", "100");
        track.setAttribute("aria-valuenow", String(percent));
        const bar = document.createElement("span");
        bar.className = "job-center-bar";
        bar.style.width = `${percent}%`;
        track.append(bar);
        button.append(track);
      }
      item.append(button);
      if (["RUNNING", "QUEUED", "WAITING_FOR_AF3_UPLOAD"].includes(job.status)) {
        const cancel = document.createElement("button");
        cancel.type = "button";
        cancel.className = "job-center-cancel";
        cancel.textContent = copy().cancel;
        cancel.setAttribute("aria-label", `${copy().cancel}: ${jobLabel(job)} (${job.id})`);
        cancel.addEventListener("click", async () => {
          cancel.disabled = true;
          try {
            const response = await fetch(`/api/jobs/${encodeURIComponent(job.id)}/cancel`, { method: "POST" });
            if (!response.ok) throw new Error(`HTTP ${response.status}`);
            await loadJobs();
          } catch {
            state.error = copy().cancelError;
            render();
          }
        });
        item.append(cancel);
      }
      list.append(item);
    });
    container.append(list);
  }

  function render() {
    const labels = copy();
    const trigger = document.getElementById("global-job-center-trigger");
    if (!trigger) return;
    const { active, recent } = partitionJobs(state.jobs);
    trigger.querySelector("[data-job-center-label]").textContent = labels.tasks;
    trigger.classList.toggle("has-active", active.length > 0);
    trigger.setAttribute("aria-label", active.length
      ? `${labels.title}: ${active.length}`
      : labels.title);
    const count = trigger.querySelector(".job-center-count");
    count.hidden = active.length === 0;
    count.textContent = String(active.length);

    const dialog = document.getElementById("global-job-center");
    dialog.querySelector("[data-job-center-title]").textContent = labels.title;
    dialog.querySelector("[data-job-center-intro]").textContent = labels.intro;
    dialog.querySelector("[data-job-center-active]").textContent = labels.active;
    dialog.querySelector("[data-job-center-recent]").textContent = labels.recent;
    const close = dialog.querySelector("[data-job-center-close]");
    close.setAttribute("aria-label", labels.close);
    close.title = labels.close;
    renderList(dialog.querySelector("[data-job-center-active-list]"), active, labels.emptyActive);
    renderList(dialog.querySelector("[data-job-center-recent-list]"), recent, labels.emptyRecent);
    const error = dialog.querySelector("[data-job-center-error]");
    error.hidden = !state.error;
    error.textContent = state.error || "";
  }

  async function loadJobs() {
    if (state.loading) return;
    if (document.visibilityState === "hidden") return;
    state.loading = true;
    try {
      const response = await fetch("/api/jobs", {
        cache: "no-store",
        headers: { Accept: "application/json" },
      });
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const jobs = await response.json();
      state.jobs = Array.isArray(jobs) ? jobs : [];
      state.error = "";
    } catch {
      state.error = copy().loadError;
    } finally {
      state.loading = false;
      render();
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", build, { once: true });
  } else {
    build();
  }
})();
