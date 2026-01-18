/**
 * Variable Name Mappings - Build Output Only
 * FULL RENAME including storage keys (no migration needed - single user)
 * 
 * NOTE: Some storage keys are NOT renamed because popup-animation.js
 * (physics file) is not processed by build and uses original names
 */

module.exports = {
    variables: [
        // ============================================
        // STORAGE KEYS - Renamed (where safe)
        // SKIPPED: moleculeCount, moleculeScale, cursorGravityStrength, 
        // initialVelocity, uiBlur, uiOpacity - these are used by physics file
        // ============================================
        ['enabled', 'jedimode', 'all'],
        ['performanceMode', 'hyperdrive', 'all'],
        ['devMode', 'sithmode', 'all'],
        ['enableAIFlagControl', 'forceguidance', 'all'],
        ['saveSizePerImage', 'savesizelightsaber', 'all'],
        ['saveSizeBySMILES', 'savesizeholocron', 'all'],
        ['useStereochemistry', 'use3dforce', 'all'],
        ['sdShowCarbons', 'showcarbonz', 'all'],
        ['sdAromaticRings', 'aromaticringz', 'all'],
        ['sdTheme', 'viztheme', 'all'],
        ['sdAutoAdapt', 'autoadaptmode', 'all'],
        ['sdRotate', 'rotatemode', 'all'],
        ['sdAverageSize', 'avgsize', 'all'],
        ['sdGradientColors', 'gradientcolors', 'all'],
        ['sdScaleByWeight', 'scalebyweight', 'all'],
        ['searchCompounds', 'findcompounds', 'all'],
        ['searchBiomolecules', 'findproteins', 'all'],
        ['searchMinerals', 'findcrystals', 'all'],
        ['proteinRemoveWhiteBg', 'removebg', 'all'],
        ['molviewBioAssembly', 'bioassembly', 'all'],
        ['bioAssemblyViewer', 'bioviewer', 'all'],
        ['bioAssemblyGraphics', 'biographics', 'all'],
        ['bioAssemblyPdbProvider', 'pdbprovider', 'all'],
        ['rendererEngine', 'renderengine', 'all'],
        // SKIPPED - used by popup-animation.js (not processed by build):
        // ['uiBlur', 'bluramt', 'all'],
        // ['uiOpacity', 'opacityamt', 'all'],
        // ['moleculeCount', 'molcount', 'all'],
        // ['moleculeScale', 'molscale', 'all'],
        // ['cursorGravityStrength', 'gravitypull', 'all'],
        // ['initialVelocity', 'startspeed', 'all'],
        ['collapsedSections', 'hiddensections', 'all'],
        ['compound3DSize', 'compound3dsize', 'all'],
        ['mineral3DSize', 'mineral3dsize', 'all'],


        // ============================================
        // CHEMISTRY THEMED - Semi-perfect (25%)
        // ============================================
        ['smilesCache', 'moleculeCache', 'content'],
        ['smilesCacheLoaded', 'moleculeCacheReady', 'content'],
        ['loadSmilesCache', 'loadMoleculeCache', 'content'],
        ['saveSmilesCache', 'saveMoleculeCache', 'content'],
        ['renderedImageCache', 'structureRenderCache', 'content'],
        ['SMILES_CACHE_KEY', 'MOLECULE_CACHE_KEY', 'content'],
        ['SMILES_CACHE_MAX_SIZE', 'MOLECULE_CACHE_LIMIT', 'content'],
        ['generateImageCacheKey', 'genStructureKey', 'content'],
        ['ELEMENT_MARKERS', 'ATOM_COLOR_MARKERS', 'content'],
        ['THEME_COLORS', 'ELEMENT_COLOR_SCHEMES', 'content'],
        ['moleculeData', 'compoundInfo', 'content'],
        ['buildMolViewEmbedUrl', 'buildMolecularViewerUrl', 'content'],
        ['buildBioImageUrl', 'buildProteinImageUrl', 'content'],
        ['wrapChemicalFormulas', 'wrapMolecularFormulas', 'content'],
        ['renderedMolecules', 'renderedCompounds', 'content'],
        ['renderClientSide', 'renderMoleculeLocally', 'content'],
        ['renderBiomolecule2D', 'renderProtein2D', 'content'],
        ['stripStereochemistry', 'removeChiralCenters', 'content'],
        ['clearAllCaches', 'executechemistryorder66', 'content'],
        ['reRenderAllMolecules', 'reactioncomplete420', 'content'],
        ['MAX_IMAGES_BEING_RENDERED', 'maxbondsinprogress', 'content'],

        // ============================================
        // TRASHY MEMES - Star Wars (40%)
        // ============================================
        ['smilesBridge', 'yodabridge', 'content'],
        ['fetchFromChemistryLaTeXServer', 'askdarthvader69', 'content'],
        ['directFetchBlob', 'usedaforce', 'content'],
        ['backgroundFetchBlob', 'padmewhenanikillsyounglings', 'content'],
        ['CHEMISTRYLATEX_CACHE_VERSION', 'deathstarversion', 'content'],
        ['cacheBustTimestamp', 'anakinisangrybecausepadmedied', 'content'],
        ['imagesBeingRendered', 'padawanzzzz', 'content'],
        ['MAX_IMAGE_CACHE_SIZE', 'maxjedicache', 'content'],
        ['getImageKey', 'getlightsaberkey', 'content'],
        ['getPageImageKey', 'ihateyouobi1', 'content'],
        ['MAX_DEFAULT_SIZE', 'youweretheoneanakin', 'content'],
        ['applyAverageSizeScaling', 'idontlikesanditscoarse', 'content'],
        ['applyScaleToImage', 'fromhereeverythingissoft', 'content'],
        ['renderSelectionAsMolecule', 'padmestilllovedanakineventhoughhedidwarcrimes', 'content'],

        // ============================================
        // TRASHY MEMES - Harry Potter (25%)
        // ============================================
        // REMOVED: getCSSFilterForTheme, getEffectiveTheme, isDarkModeEnabled, pageIsDark - all auto-adapt code deleted

        // ============================================
        // TRASHY - The Office
        // ============================================
        ['show3DViewerInline', 'dwightshows3d', 'content'],
        ['addHoverControls', 'dundermifflincontrols', 'content'],
        ['wrapImageWithSizeControls', 'thatswhatsheresized', 'content'],
        ['wrapImageWithRotationContainer', 'michaelscottspin', 'content'],
        ['reloadAllImages', 'kevindroppedchili', 'content'],
        ['reRenderSingleMolecule', 'tobynobodylikeshim', 'content'],
        ['processImageIfNeeded', 'processliketoby', 'content'],

        // ============================================
        // LAZY LOADING - Keep professional
        // ============================================
        ['setupLazyLoading', 'initializeLazyLoader', 'content'],
        ['lazyReRenderMolecules', 'deferredRenderUpdate', 'content'],
        ['lazyReRenderObserver', 'intersectionObserverInstance', 'content'],
        ['cleanupObservers', 'disconnectAllObservers', 'content'],
        ['triggerVisibleImagesReRender', 'refreshVisibleElements', 'content'],

        // ============================================
        // TRASHY - Minecraft
        // ============================================
        ['BASE_SIZE', 'cobblestonebase', 'content'],
        ['DEFAULT_WIDTH', 'stevewidth', 'content'],
        ['DEFAULT_HEIGHT', 'steveheight', 'content'],
        ['MAX_SIZE', 'enderdragonsiz', 'content'],
        ['MIN_SIZE', 'babyzombie', 'content'],
        ['SIZE_STEP', 'blockstep', 'content'],
        ['SIZE_SCALE_FACTOR', 'stevescale', 'content'],
        ['loadImageSize', 'loadblockdimension', 'content'],
        ['saveImageSize', 'saveblockdimension', 'content'],

        // ============================================
        // More Star Wars memes
        // ============================================
        ['pendingReRenderSettings', 'itsatrappp', 'content'],
        ['isUpdatingImages', 'isrefreshn', 'content'],
        ['applyCssOnlySettings', 'applydriponly', 'content'],
        ['addLoadingIndicator', 'showspinner2', 'content'],
        ['removeLoadingIndicator', 'hidespinner2', 'content'],
        ['pendingUpdateIndicator', 'dikiloveyou_iknow', 'content'],
        ['injectStyles', 'injectdadrip', 'content'],
        ['getOrCreatePendingSearch', 'dikiloveyou2', 'content'],
        ['pendingSearches', 'searchinprogress', 'content'],
        ['MAX_PENDING_SEARCHES', 'neversaytometheodds', 'content'],
        ['getAffectedImageTypes', 'getaffectedtypez', 'content'],
        ['log', 'bobthelogger', 'content'],
        ['LOG_PREFIX', 'bobtag', 'content'],
        ['logHistory', 'bobsnotes', 'content'],
        ['extensionContextInvalid', 'theempiredidnothingwrong', 'content'],
        ['isExtensionContextValid', 'rebelscumisalive', 'content'],
        ['activeLoads', 'wateronmyneck', 'content'],  // tracking the pressure/load on the system

        ['viewerUrl', 'viewersrc', 'content'],
        ['viewer3DIframe', 'iframe3d', 'content'],
        ['bioToggleBtn', 'biotoggle', 'content'],
        ['cubeBtn', 'thisistehway', 'content'],
        ['buttonGroup', 'btngroup', 'content'],
        ['hoverControls', 'floatybois', 'content'],

        // ============================================
        // SCAN functions
        // ============================================
        ['scanAndRender', 'senddaprobes', 'content'],
        ['scanAndRenderImmediate', 'probesnow', 'content'],
        ['scanAndRenderTimeout', 'probedelay', 'content'],
        ['lastScanTime', 'lastprobetime', 'content'],
        ['MIN_SCAN_INTERVAL', 'probecooldown', 'content'],
        ['initializeRenderer', 'bootxwing', 'content'],
        ['setupLazyReRenderObserver', 'initneowatchr', 'content'],
    ],

    // String replacements - only safe ones
    strings: [
        ['SmilesDrawer', 'MoleculeEngine', 'all'],
        ['smiles-drawer', 'molecule-engine', 'all'],
        ['smilesDrawer', 'moleculeEngine', 'all'],
    ],
};
