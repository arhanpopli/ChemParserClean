/**
 * Integrated Search Module - NO LOCAL DATABASES NEEDED!
 * 
 * This module queries ALL three major chemistry databases directly:
 * 1. RCSB PDB API - for biomolecules/proteins (250k+ entries)
 * 2. COD API - for minerals/crystals (500k+ entries)  
 * 3. PubChem API - for organic compounds (100M+ entries)
 * 
 * Then uses similar_text algorithm to find the best match across all results.
 * Returns pdbid for biomolecules, codid for minerals, smiles for compounds.
 */

(function () {
    'use strict';

    // ============================================
    // SESSION CACHE - Deduplicates API calls per page session
    // ============================================

    // In-memory cache for search results (cleared on page reload)
    const sessionCache = new Map();

    // In-memory map of pending searches (for deduplication)
    const pendingSearches = new Map();

    /**
     * Get cached result or pending search promise
     * @param {string} query - Search query
     * @returns {object|null} - Cached result or null
     */
    function getCachedResult(query) {
        const key = query.toLowerCase().trim();
        if (sessionCache.has(key)) {
            console.log(`[IntegratedSearch] 🎯 CACHE HIT: "${query}"`);
            return sessionCache.get(key);
        }
        return null;
    }

    /**
     * Store result in session cache
     * @param {string} query - Search query
     * @param {object} result - Search result
     */
    function setCachedResult(query, result) {
        const key = query.toLowerCase().trim();
        sessionCache.set(key, result);
        console.log(`[IntegratedSearch] 📦 Cached result for: "${query}"`);
    }

    // ============================================
    // UTILITY FUNCTIONS
    // ============================================

    function ucfirst(str) {
        if (str === undefined) return str;
        return str.charAt(0).toUpperCase() + str.slice(1);
    }

    function humanize(str) {
        if (str === undefined) return str;
        return str.replace(/(\b[A-Z]+\b)/g, function (word) {
            return word.toLowerCase();
        });
    }

    // PHP similar_text algorithm for fuzzy matching
    function similar_text(first, second, percent) {
        if (first === null || second === null || typeof first === 'undefined' || typeof second === 'undefined') {
            return 0;
        }

        first += '';
        second += '';

        var pos1 = 0, pos2 = 0, max = 0;
        var firstLength = first.length, secondLength = second.length;
        var p, q, l, sum;

        for (p = 0; p < firstLength; p++) {
            for (q = 0; q < secondLength; q++) {
                for (l = 0; (p + l < firstLength) && (q + l < secondLength) && (first.charAt(p + l) === second.charAt(q + l)); l++);
                if (l > max) {
                    max = l;
                    pos1 = p;
                    pos2 = q;
                }
            }
        }

        sum = max;

        if (sum) {
            if (pos1 && pos2) {
                sum += similar_text(first.substr(0, pos1), second.substr(0, pos2));
            }
            if ((pos1 + max < firstLength) && (pos2 + max < secondLength)) {
                sum += similar_text(first.substr(pos1 + max, firstLength - pos1 - max),
                    second.substr(pos2 + max, secondLength - pos2 - max));
            }
        }

        return !percent ? sum : (sum * 200) / (firstLength + secondLength);
    }

    // ============================================
    // AUTOCOMPLETE BUILDER
    // ============================================

    class AutocompleteBuilder {
        constructor(array, key) {
            this.array = array || [];
            this.key = key;
        }

        sort(str, minsim, length) {
            str = str.toLowerCase();
            var cpy = this.array.slice(0);
            for (var i = 0; i < cpy.length; i++) {
                cpy[i].similarity = similar_text(str, cpy[i][this.key].toLowerCase(), true);

                // Add 100 to similarity if first characters match
                if (cpy[i][this.key].toLowerCase().indexOf(str) === 0) {
                    cpy[i].similarity += 100;
                }

                // CRITICAL: Add 200 if the result CONTAINS the exact query word
                // This ensures "rhinovirus-B5" beats "Novirus" when searching "rhinovirus"
                if (cpy[i][this.key].toLowerCase().includes(str)) {
                    cpy[i].similarity += 200;
                }

                if (minsim && cpy[i].similarity < minsim) {
                    cpy.splice(i, 1);
                    i--;
                }
            }

            return cpy.sort((a, b) => b.similarity - a.similarity)
                .slice(0, length ? length : cpy.length);
        }
    }

    // ============================================
    // FETCH HELPER (with CSP bypass)
    // ============================================

    async function fetchUrl(url, timeout = 30000) {
        if (typeof chrome !== 'undefined' && chrome.runtime && chrome.runtime.sendMessage) {
            try {
                return await new Promise((resolve, reject) => {
                    const timeoutId = setTimeout(() => reject(new Error('Timeout')), timeout);
                    chrome.runtime.sendMessage(
                        { type: 'FETCH_TEXT', url: url, options: {} },
                        (response) => {
                            clearTimeout(timeoutId);
                            if (chrome.runtime.lastError) {
                                directFetch(url, timeout).then(resolve).catch(reject);
                                return;
                            }
                            if (response && response.success) {
                                resolve(response.data);
                            } else {
                                reject(new Error(response?.error || 'Background fetch failed'));
                            }
                        }
                    );
                });
            } catch (e) {
                return directFetch(url, timeout);
            }
        }
        return directFetch(url, timeout);
    }

    async function directFetch(url, timeout = 30000) {
        const controller = new AbortController();
        const timeoutId = setTimeout(() => controller.abort(), timeout);

        try {
            const response = await fetch(url, { signal: controller.signal });
            clearTimeout(timeoutId);
            if (!response.ok) throw new Error(`HTTP ${response.status}`);
            return await response.text();
        } catch (error) {
            clearTimeout(timeoutId);
            throw error;
        }
    }

    // ============================================
    // CIF PARSER
    // ============================================

    function parseCIFData(cifContent) {
        const data = {
            chemical_name: null,
            formula: null,
            canonical_smiles: null,
            isomeric_smiles: null
        };

        if (!cifContent) return data;

        const nameMatch = cifContent.match(/_chemical_name_common\s+['"]?([^'"\r\n]+)['"]?/i) ||
            cifContent.match(/_chemical_name_mineral\s+['"]?([^'"\r\n]+)['"]?/i);
        if (nameMatch) data.chemical_name = nameMatch[1].trim();

        const formulaMatch = cifContent.match(/_chemical_formula_sum\s+['"]?([^'"\r\n]+)['"]?/i);
        if (formulaMatch) data.formula = formulaMatch[1].trim().replace(/['"]/g, '');

        const canonicalMatch = cifContent.match(/_chemical_smiles_canonical\s+['"]?([^'"\r\n]+)['"]?/i);
        if (canonicalMatch) data.canonical_smiles = canonicalMatch[1].trim();

        const isomericMatch = cifContent.match(/_chemical_smiles_isomeric\s+['"]?([^'"\r\n]+)['"]?/i);
        if (isomericMatch) data.isomeric_smiles = isomericMatch[1].trim();

        return data;
    }

    // ============================================
    // MAIN SEARCH FUNCTION - Queries ALL APIs!
    // ============================================

    async function integratedSearch(query, options = {}) {
        console.log(`[IntegratedSearch] 🔍 Searching ALL databases for: "${query}"`);

        if (!query || query.trim() === '') {
            return { error: 'No query provided' };
        }

        query = query.trim();
        const lowerQuery = query.toLowerCase();
        const cacheKey = lowerQuery;

        // ============ CHECK SESSION CACHE FIRST ============
        const cachedResult = getCachedResult(query);
        if (cachedResult) {
            return cachedResult;
        }

        // ============ CHECK IF SEARCH IS ALREADY IN PROGRESS ============
        if (pendingSearches.has(cacheKey)) {
            console.log(`[IntegratedSearch] 🔄 Waiting for pending search: "${query}"`);
            return pendingSearches.get(cacheKey);
        }

        // Create a promise for this search
        const searchPromise = (async () => {
            try {
                const result = await doIntegratedSearch(query, lowerQuery, options);
                // Cache successful results
                if (result && !result.error) {
                    setCachedResult(query, result);
                }
                return result;
            } finally {
                pendingSearches.delete(cacheKey);
            }
        })();

        pendingSearches.set(cacheKey, searchPromise);
        return searchPromise;
    }

    // Actual search implementation (called by integratedSearch)
    async function doIntegratedSearch(query, lowerQuery, options = {}) {
        // ============ 0. DIRECT ID DETECTION (pdbid / codid / cid) ============
        const directId = (() => {
            // PDB IDs: 4 chars, starting with a digit (canonical) or prefixed
            const pdbMatch = query.match(/^pdb[:\s-]?([0-9][A-Za-z0-9]{3})$/i) || query.match(/^([0-9][A-Za-z0-9]{3})$/);
            if (pdbMatch) return { type: 'pdb', value: pdbMatch[1].toUpperCase() };

            // COD IDs: numeric, usually 4-7 digits, or prefixed
            const codMatch = query.match(/^cod[:\s-]?(\d{4,7})$/i);
            if (codMatch) return { type: 'cod', value: codMatch[1] };

            // PubChem CID: numeric, or prefixed
            const cidMatch = query.match(/^cid[:\s-]?(\d+)$/i);
            if (cidMatch) return { type: 'cid', value: cidMatch[1] };

            // Bare numbers: prefer COD for long IDs, otherwise CID
            if (/^\d+$/.test(query)) {
                if (query.length >= 7) return { type: 'cod', value: query };
                return { type: 'cid', value: query };
            }
            return null;
        })();

        if (directId) {
            console.log(`[IntegratedSearch] 🎯 Direct ID detected: ${directId.type.toUpperCase()} ${directId.value}`);
            const result = {
                query,
                corrected_query: query,
                name: query,
                chemical_name: null,
                formula: null,
                canonical_smiles: null,
                isomeric_smiles: null,
                sdf: null,
                source_url: null,
                primary_type: 'compound',
                sub_type: null,
                embed_url: null
            };

            try {
                if (directId.type === 'pdb') {
                    const pdbId = directId.value;
                    result.primary_type = 'biomolecule';
                    result.pdbid = pdbId;
                    result.source_url = `https://www.rcsb.org/structure/${pdbId}`;
                    result.pdb_url = `https://files.rcsb.org/view/${pdbId}.pdb`;
                    result.embed_url = `https://embed.molview.org/v1/?pdbid=${pdbId}`;
                    result.image_url = `https://cdn.rcsb.org/images/structures/${pdbId.toLowerCase()}_model-1.jpeg`;
                    console.log(`[IntegratedSearch] ✅ Direct PDB: ${pdbId}`);
                    return result;
                }

                if (directId.type === 'cod') {
                    const codid = directId.value;
                    result.primary_type = 'mineral';
                    result.codid = codid;
                    result.source_url = `http://www.crystallography.net/cod/${codid}.html`;
                    result.cif_url = `https://www.crystallography.net/cod/${codid}.cif`;
                    result.embed_url = `https://embed.molview.org/v1/?codid=${codid}`;

                    try {
                        const cifContent = await fetchUrl(result.cif_url);
                        if (!cifContent.includes('404')) {
                            const cifData = parseCIFData(cifContent);
                            result.chemical_name = cifData.chemical_name;
                            result.formula = cifData.formula;
                            result.canonical_smiles = cifData.canonical_smiles;
                            result.isomeric_smiles = cifData.isomeric_smiles;
                        }
                    } catch (e) {
                        console.log(`[IntegratedSearch] CIF fetch failed for COD ${codid}: ${e.message}`);
                    }
                    console.log(`[IntegratedSearch] ✅ Direct COD: ${codid}`);
                    return result;
                }

                if (directId.type === 'cid') {
                    const cid = directId.value;
                    result.primary_type = 'compound';
                    result.cid = cid;
                    result.embed_url = `https://embed.molview.org/v1/?cid=${cid}`;
                    result.source_url = `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`;
                    try {
                        const propsUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/Title,SMILES,ConnectivitySMILES,MolecularFormula,IUPACName/JSON`;
                        const propsData = JSON.parse(await fetchUrl(propsUrl));
                        if (propsData.PropertyTable && propsData.PropertyTable.Properties) {
                            const props = propsData.PropertyTable.Properties[0];
                            result.name = props.Title || result.name;
                            result.chemical_name = props.IUPACName || props.Title;
                            result.formula = props.MolecularFormula;
                            result.isomeric_smiles = props.SMILES;
                            result.canonical_smiles = props.ConnectivitySMILES || props.SMILES;
                        }
                    } catch (e) {
                        console.log(`[IntegratedSearch] PubChem property fetch failed for CID ${cid}: ${e.message}`);
                    }
                    console.log(`[IntegratedSearch] ✅ Direct CID: ${cid}`);
                    return result;
                }
            } catch (e) {
                console.warn(`[IntegratedSearch] Direct ID handling failed: ${e.message}`);
            }
        }

        let mix = [];

        // ============ 1. RCSB PDB API - Search for proteins/biomolecules ============
        console.log('[IntegratedSearch] 🧬 Querying RCSB PDB API...');
        try {
            // Use RCSB text search API
            const rcsbSearchUrl = `https://search.rcsb.org/rcsbsearch/v2/query`;
            const rcsbQuery = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": query
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "return_all_hits": false,
                    "paginate": {
                        "start": 0,
                        "rows": 10
                    }
                }
            };

            const rcsbResponse = await fetch(rcsbSearchUrl, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(rcsbQuery)
            });

            if (rcsbResponse.ok) {
                const rcsbData = await rcsbResponse.json();
                if (rcsbData.result_set && rcsbData.result_set.length > 0) {
                    // Get details for top results
                    for (let i = 0; i < Math.min(5, rcsbData.result_set.length); i++) {
                        const pdbid = rcsbData.result_set[i].identifier;
                        try {
                            const detailsUrl = `https://data.rcsb.org/rest/v1/core/entry/${pdbid}`;
                            const detailsResp = await fetch(detailsUrl);
                            if (detailsResp.ok) {
                                const details = await detailsResp.json();
                                const title = details.struct?.title || details.rcsb_entry_info?.title || pdbid;
                                mix.push({
                                    name: title,
                                    label: title,
                                    pdbids: [pdbid],
                                    source_type: 'rcsb',
                                    score: rcsbData.result_set[i].score || 0
                                });
                                console.log(`[IntegratedSearch] Found RCSB: ${title} (${pdbid})`);
                            }
                        } catch (e) {
                            console.log(`[IntegratedSearch] Failed to get details for ${pdbid}`);
                        }
                    }
                }
            }
        } catch (e) {
            console.log(`[IntegratedSearch] RCSB search failed: ${e.message}`);
        }

        // ============ 2. COD API - Search for minerals/crystals ============
        console.log('[IntegratedSearch] 💎 Querying COD API...');
        try {
            const codSearchUrl = `https://www.crystallography.net/cod/result?text=${encodeURIComponent(query)}&format=json`;
            const codDataStr = await fetchUrl(codSearchUrl);
            const codData = JSON.parse(codDataStr);

            if (codData && Array.isArray(codData) && codData.length > 0) {
                // Take top 10 results
                for (let i = 0; i < Math.min(10, codData.length); i++) {
                    const entry = codData[i];
                    const codid = entry.file || entry.id;
                    const name = entry.commonname || entry.chemname || entry.formula || `COD ${codid}`;

                    mix.push({
                        name: name,
                        label: ucfirst(name),
                        codid: codid,
                        source_type: 'cod'
                    });
                    console.log(`[IntegratedSearch] Found COD: ${name} (${codid})`);
                }
            }
        } catch (e) {
            console.log(`[IntegratedSearch] COD search failed: ${e.message}`);
        }

        // ============ 3. PUBCHEM - Direct lookup + Autocomplete ============
        console.log('[IntegratedSearch] 🧪 Querying PubChem API...');

        // 3a. Direct lookup (exact match)
        try {
            const directUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(query)}/cids/JSON`;
            const directDataStr = await fetchUrl(directUrl);
            if (!directDataStr.includes('PUGREST.NotFound')) {
                const directData = JSON.parse(directDataStr);
                if (directData.IdentifierList && directData.IdentifierList.CID && directData.IdentifierList.CID.length > 0) {
                    console.log(`[IntegratedSearch] ✅ Direct PubChem hit: CID ${directData.IdentifierList.CID[0]}`);
                    mix.push({
                        name: query,
                        label: ucfirst(humanize(query)),
                        PubChem_name: query,
                        source_type: 'pubchem_direct',
                        similarity: 300  // High priority for exact matches
                    });
                }
            }
        } catch (e) {
            console.log(`[IntegratedSearch] PubChem direct lookup: not found`);
        }

        // 3b. Autocomplete
        try {
            const autocpUrl = `https://pubchem.ncbi.nlm.nih.gov/pcautocp/pcautocp.cgi?dict=pc_compoundnames&n=10&q=${encodeURIComponent(query)}`;
            const autocpDataStr = await fetchUrl(autocpUrl);
            const autocpData = JSON.parse(autocpDataStr);

            if (autocpData.autocp_array) {
                for (let name of autocpData.autocp_array) {
                    mix.push({
                        name: name,
                        label: ucfirst(humanize(name)),
                        PubChem_name: name,
                        source_type: 'pubchem'
                    });
                }
                console.log(`[IntegratedSearch] PubChem autocomplete: ${autocpData.autocp_array.length} results`);
            }
        } catch (e) {
            console.log(`[IntegratedSearch] PubChem autocomplete failed: ${e.message}`);
        }

        // ============ 4. SORT ALL RESULTS BY SIMILARITY ============
        console.log(`[IntegratedSearch] Total results from all APIs: ${mix.length}`);

        if (mix.length === 0) {
            return { error: 'No results found in any database', query: query };
        }

        const finalSorter = new AutocompleteBuilder(mix, "label");
        let sortedMix = finalSorter.sort(query);

        // Remove duplicates
        for (var i = 0; i < sortedMix.length; i++) {
            if (i < sortedMix.length - 1) {
                if (sortedMix[i].label && sortedMix[i + 1].label &&
                    sortedMix[i].label.toLowerCase() === sortedMix[i + 1].label.toLowerCase()) {
                    sortedMix[i].codid = sortedMix[i].codid || sortedMix[i + 1].codid;
                    sortedMix[i].pdbids = sortedMix[i].pdbids || sortedMix[i + 1].pdbids;
                    sortedMix[i].PubChem_name = sortedMix[i].PubChem_name || sortedMix[i + 1].PubChem_name;
                    sortedMix.splice(i + 1, 1);
                }
            }
        }

        // ============ 5. FIND BEST MATCH ============
        let bestMatch = sortedMix[0];
        const queryLower = query.toLowerCase();

        // Prefer exact matches
        for (let i = 0; i < sortedMix.length; i++) {
            if (sortedMix[i].name && sortedMix[i].name.toLowerCase() === queryLower) {
                bestMatch = sortedMix[i];
                console.log(`[IntegratedSearch] ✅ Exact name match: ${bestMatch.name}`);
                break;
            }
            if (sortedMix[i].label && sortedMix[i].label.toLowerCase() === queryLower) {
                bestMatch = sortedMix[i];
                console.log(`[IntegratedSearch] ✅ Exact label match: ${bestMatch.label}`);
                break;
            }
        }

        console.log(`[IntegratedSearch] 🏆 Best match: ${bestMatch.label || bestMatch.name} (source: ${bestMatch.source_type}, pdbid: ${bestMatch.pdbids ? bestMatch.pdbids[0] : 'none'}, codid: ${bestMatch.codid || 'none'})`);

        // ============ 6. BUILD RESULT ============
        const result = {
            query: query,
            corrected_query: bestMatch.label,
            name: bestMatch.name || bestMatch.label,
            chemical_name: null,
            formula: null,
            canonical_smiles: null,
            isomeric_smiles: null,
            sdf: null,
            source_url: null,
            primary_type: 'compound',
            sub_type: null,
            embed_url: null
        };

        // ============ BIOMOLECULE (has pdbids) ============
        if (bestMatch.pdbids && bestMatch.pdbids.length > 0) {
            const pdbId = bestMatch.pdbids[0];
            result.primary_type = 'biomolecule';
            result.pdbid = pdbId;
            result.source_url = `https://www.rcsb.org/structure/${pdbId}`;
            result.pdb_url = `https://files.rcsb.org/view/${pdbId}.pdb`;
            result.embed_url = `https://embed.molview.org/v1/?pdbid=${pdbId}`;
            result.image_url = `https://cdn.rcsb.org/images/structures/${pdbId.toLowerCase()}_model-1.jpeg`;

            console.log(`[IntegratedSearch] ✅ BIOMOLECULE: ${bestMatch.name} (PDB: ${pdbId})`);

            // ============ MINERAL (has codid) ============
        } else if (bestMatch.codid) {
            const codid = bestMatch.codid;
            result.primary_type = 'mineral';
            result.codid = codid;
            result.source_url = `http://www.crystallography.net/cod/${codid}.html`;
            result.cif_url = `https://www.crystallography.net/cod/${codid}.cif`;
            result.embed_url = `https://embed.molview.org/v1/?codid=${codid}`;

            console.log(`[IntegratedSearch] ✅ MINERAL: ${bestMatch.name} (COD: ${codid})`);

            // Fetch CIF for formula and SMILES
            try {
                const cifContent = await fetchUrl(result.cif_url);
                if (!cifContent.includes('404')) {
                    const cifData = parseCIFData(cifContent);
                    result.chemical_name = cifData.chemical_name;
                    result.formula = cifData.formula;
                    result.canonical_smiles = cifData.canonical_smiles;
                    result.isomeric_smiles = cifData.isomeric_smiles;
                    console.log(`[IntegratedSearch] CIF data: formula=${result.formula}, smiles=${result.canonical_smiles || 'none'}`);
                }
            } catch (e) {
                console.log(`[IntegratedSearch] CIF fetch failed: ${e.message}`);
            }

            // Try PubChem for SMILES if not in CIF
            if (!result.canonical_smiles) {
                const mineralToChemical = {
                    'Quartz': 'silicon dioxide', 'Calcite': 'calcium carbonate',
                    'Gypsum': 'calcium sulfate', 'Halite': 'sodium chloride',
                    'Fluorite': 'calcium fluoride', 'Diamond': 'carbon', 'Graphite': 'carbon'
                };
                let searchName = mineralToChemical[bestMatch.name] || bestMatch.name;

                try {
                    const cidUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(searchName)}/cids/JSON`;
                    const cidDataStr = await fetchUrl(cidUrl);
                    if (!cidDataStr.includes('PUGREST.NotFound')) {
                        const cidData = JSON.parse(cidDataStr);
                        if (cidData.IdentifierList && cidData.IdentifierList.CID) {
                            const cid = cidData.IdentifierList.CID[0];
                            const propsUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/SMILES,ConnectivitySMILES,MolecularFormula/JSON`;
                            const propsData = JSON.parse(await fetchUrl(propsUrl));
                            if (propsData.PropertyTable && propsData.PropertyTable.Properties) {
                                const props = propsData.PropertyTable.Properties[0];
                                result.isomeric_smiles = props.SMILES;
                                result.canonical_smiles = props.ConnectivitySMILES || props.SMILES;
                                result.formula = result.formula || props.MolecularFormula;
                            }
                        }
                    }
                } catch (e) { }
            }

            // ============ COMPOUND (PubChem) ============
        } else {
            result.primary_type = 'compound';
            let cid = null;

            // Get CID
            try {
                const cidUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(bestMatch.name || bestMatch.PubChem_name)}/cids/JSON`;
                const cidDataStr = await fetchUrl(cidUrl);
                if (!cidDataStr.includes('PUGREST.NotFound')) {
                    const cidData = JSON.parse(cidDataStr);
                    if (cidData.IdentifierList && cidData.IdentifierList.CID) {
                        cid = cidData.IdentifierList.CID[0];
                    }
                }
            } catch (e) { }

            // Try original query if failed
            if (!cid) {
                try {
                    const originalCidUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(query)}/cids/JSON`;
                    const originalCidDataStr = await fetchUrl(originalCidUrl);
                    if (!originalCidDataStr.includes('PUGREST.NotFound')) {
                        const originalCidData = JSON.parse(originalCidDataStr);
                        if (originalCidData.IdentifierList && originalCidData.IdentifierList.CID) {
                            cid = originalCidData.IdentifierList.CID[0];
                        }
                    }
                } catch (e) { }
            }

            if (cid) {
                result.cid = cid;
                result.embed_url = `https://embed.molview.org/v1/?cid=${cid}`;
                result.source_url = `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`;

                // Fetch properties
                try {
                    const propsUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/Title,SMILES,ConnectivitySMILES,MolecularFormula,IUPACName/JSON`;
                    const propsData = JSON.parse(await fetchUrl(propsUrl));
                    if (propsData.PropertyTable && propsData.PropertyTable.Properties) {
                        const props = propsData.PropertyTable.Properties[0];
                        result.name = props.Title || result.name;
                        result.chemical_name = props.IUPACName || props.Title;
                        result.formula = props.MolecularFormula;
                        result.isomeric_smiles = props.SMILES;
                        result.canonical_smiles = props.ConnectivitySMILES || props.SMILES;
                        console.log(`[IntegratedSearch] ✅ COMPOUND: ${result.name} (CID: ${cid})`);
                    }
                } catch (e) { }
            }
        }

        console.log(`[IntegratedSearch] ✅ Complete: {type: '${result.primary_type}', name: '${result.name}', pdbid: '${result.pdbid || 'N/A'}', codid: '${result.codid || 'N/A'}', smiles: '${result.canonical_smiles ? result.canonical_smiles.substring(0, 20) + '...' : 'N/A'}'}`);
        return result;
    }

    // ============================================
    // EXPORT
    // ============================================

    window.IntegratedSearch = {
        search: integratedSearch,
        clearCache: function () {
            sessionCache.clear();
            pendingSearches.clear();
            console.log('[IntegratedSearch] 🗑️ Cache cleared!');
        }
    };

    console.log('[IntegratedSearch] ✅ Module loaded - Queries RCSB, COD, and PubChem APIs directly. NO local databases needed!');

})();
