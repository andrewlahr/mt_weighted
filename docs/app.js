/* ==========================================================================
   Flow–Recruitment FPCA Results Viewer — site loader
   Reads sites_index.json, populates dropdown, swaps images per site.
   ========================================================================== */

(function () {
  'use strict';

  const select   = document.getElementById('site-select');
  const loading  = document.getElementById('loading');
  const results  = document.getElementById('results');
  const meta     = document.getElementById('site-meta');
  const metaY    = document.getElementById('meta-years');
  const metaF    = document.getElementById('meta-fpcs');
  const metaR    = document.getElementById('meta-r2');
  const imgPanel = document.getElementById('img-panels');
  const imgSumm  = document.getElementById('img-summary');
  const imgFpcs  = document.getElementById('img-allfpcs');
  const imgYrFpc = document.getElementById('img-yearfpc');

  // ----- Persist last selected site in URL hash and localStorage --------
  const STORAGE_KEY = 'fpca-last-site';

  function getInitialSite(siteIds) {
    const hash = window.location.hash.replace('#', '');
    if (hash && siteIds.includes(hash)) return hash;
    const stored = localStorage.getItem(STORAGE_KEY);
    if (stored && siteIds.includes(stored)) return stored;
    return siteIds[0] || null;
  }

  function setSelected(siteId) {
    if (!siteId) return;
    select.value = siteId;
    window.history.replaceState(null, '', '#' + siteId);
    localStorage.setItem(STORAGE_KEY, siteId);
    loadSite(siteId);
  }

  // ----- Load and display a single site ---------------------------------
  async function loadSite(siteId) {
    loading.textContent = 'Loading ' + siteId + '…';
    loading.style.display = 'block';
    results.style.display = 'none';
    meta.style.display = 'none';

    try {
      const metaResp = await fetch('assets/sites/' + siteId + '/metadata.json',
                                   { cache: 'no-cache' });
      if (!metaResp.ok) throw new Error('metadata.json not found for ' + siteId);
      const m = await metaResp.json();

      // Update meta bar
      metaY.textContent = m.n_years ?? '—';
      metaF.textContent = m.n_fpc_retained ?? '—';
      metaR.textContent = (m.r_squared != null) ? m.r_squared.toFixed(3) : '—';
      meta.style.display = 'flex';

      // Swap images (cache-bust with timestamp on first load only)
      const base = 'assets/sites/' + siteId + '/';
      imgPanel.src = base + 'page1_panels.png';
      imgSumm.src  = base + 'page2_summary.png';
      imgFpcs.src  = base + 'page3_all_fpcs.png';
      imgYrFpc.src = base + 'page4_year_fpc_heatmap.png';

      // Update doc title
      document.title = m.site + ' — Flow–Recruitment FPCA';

      loading.style.display = 'none';
      results.style.display = 'block';
    } catch (err) {
      console.error(err);
      loading.innerHTML =
        '<div class="error-box">Could not load site <strong>' + siteId +
        '</strong>. Check that <code>assets/sites/' + siteId +
        '/metadata.json</code> exists.</div>';
    }
  }

  // ----- Initialise on page load ----------------------------------------
  async function init() {
    try {
      const resp = await fetch('sites_index.json', { cache: 'no-cache' });
      if (!resp.ok) throw new Error('sites_index.json not found');
      const idx = await resp.json();
      const sites = idx.sites || [];

      if (sites.length === 0) {
        select.innerHTML = '<option value="">No sites available</option>';
        loading.innerHTML =
          '<div class="error-box">No sites found in <code>sites_index.json</code>. ' +
          'Run <code>WorkingFPCA.R</code> to generate site results.</div>';
        return;
      }

      // Populate dropdown
      select.innerHTML = sites.map(s =>
        '<option value="' + s.id + '">' + s.site + '</option>'
      ).join('');

      // Bind change handler
      select.addEventListener('change', () => setSelected(select.value));

      // Bind hash change handler (for direct links)
      window.addEventListener('hashchange', () => {
        const h = window.location.hash.replace('#', '');
        if (h && h !== select.value) setSelected(h);
      });

      // Load the initial site
      const ids = sites.map(s => s.id);
      const initial = getInitialSite(ids);
      if (initial) setSelected(initial);
    } catch (err) {
      console.error(err);
      select.innerHTML = '<option value="">(error)</option>';
      loading.innerHTML =
        '<div class="error-box">Could not load <code>sites_index.json</code>. ' +
        'Make sure <code>WorkingFPCA.R</code> has been run and the file exists ' +
        'in the <code>docs/</code> directory.</div>';
    }
  }

  init();
})();
