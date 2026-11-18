class Chemistry {

    private static final Map ALIASES = [
        '10x'     : '10X',
        '10xv2'   : '10Xv2',
        '10xv3'   : '10Xv3',
        'droplet' : 'Droplet',
        'scalebio': 'ScaleBio',
        'smart'   : 'Smart',
        'smartseq': 'Smart',
        'pipseq'  : 'PIPseq'
    ]

    private static final Map LEGACY = [
        '10X'     : 'star_10X',
        '10Xv2'   : 'star_10X',
        '10Xv3'   : 'star_10X',
        'Droplet' : 'star_droplet',
        'ScaleBio': 'star_scalebio',
        'Smart'   : 'star_smart',
        'PIPseq'  : 'star_pipseq'
    ]

    private static final Map GEOMETRY_OPTIONS = [
        'soloCBstart' : 'cbStart',
        'soloCBlen'   : 'cbLen',
        'soloUMIstart': 'umiStart',
        'soloUMIlen'  : 'umiLen'
    ]

    static Map resolve(Object name, Map params) {
        def requested = name?.toString()?.trim()
        if (requested == null || requested.isEmpty()) {
            throw new IllegalArgumentException("chemistry is not set")
        }
        def definitions = params.chemistries
        if (!(definitions instanceof Map) || definitions.isEmpty()) {
            throw new IllegalArgumentException("params.chemistries is not defined")
        }
        def key = canonicalKey(requested, definitions)
        def definition = definitions[key]
        def mode = definition.mode?.toString() ?: 'CB_UMI_Simple'
        def geometry = geometryOf(key, definition, params)
        def options = optionsOf(key, definitions)
        syncGeometry(options, geometry)
        def override = legacyOverride(key, params, geometry)
        if (override != null) {
            def parsed = parseOptions(override)
            geometry = mergeGeometry(key, geometry, parsed)
            options = parsed
        }
        return [name: key, mode: mode, geometry: geometry, starOptions: options]
    }

    static String starArgs(Map resolved) {
        def options = resolved.starOptions
        if (options == null || options.isEmpty()) {
            throw new IllegalArgumentException("chemistry '${resolved.name}' has no STAR options")
        }
        return options.collect { option, value ->
            (value == null) ? "--${option}" : "--${option} ${value}"
        }.join(' ')
    }

    static String extractionArgs(Map resolved) {
        def geometry = resolved.geometry
        if (geometry == null) {
            throw new IllegalArgumentException("chemistry '${resolved.name}' has no barcode/UMI geometry, CaTCH barcode counting is not supported for it")
        }
        return "--bcStart ${geometry.cbStart - 1} --bcLength ${geometry.cbLen} --umiStart ${geometry.umiStart - 1} --umiLength ${geometry.umiLen}"
    }

    static Map parseOptions(Object text) {
        def source = text?.toString()?.trim()
        if (source == null || source.isEmpty()) {
            throw new IllegalArgumentException("STAR options must not be empty")
        }
        def options = [:]
        String current = null
        source.split(/\s+/).each { token ->
            if (token.startsWith('--')) {
                current = token.substring(2)
                if (current.isEmpty()) {
                    throw new IllegalArgumentException("malformed STAR option '${token}'")
                }
                options[current] = []
            } else if (current == null) {
                throw new IllegalArgumentException("STAR options must start with an option name, got '${token}'")
            } else {
                options[current] << token
            }
        }
        return options.collectEntries { option, values ->
            [option, values.isEmpty() ? null : values.join(' ')]
        }
    }

    private static String canonicalKey(String requested, Map definitions) {
        def direct = definitions.keySet().find { it.toString().equalsIgnoreCase(requested) }
        if (direct != null) {
            return direct.toString()
        }
        def alias = ALIASES[requested.toLowerCase()]
        if (alias != null && definitions.containsKey(alias)) {
            return alias
        }
        if (requested.toUpperCase().contains('10X') && definitions.containsKey('10X')) {
            return '10X'
        }
        throw new IllegalArgumentException("unknown chemistry '${requested}', supported chemistries are ${definitions.keySet().join(', ')}")
    }

    private static Map geometryOf(String key, Map definition, Map params) {
        def geometry = definition.geometry
        if (geometry == null) {
            return null
        }
        def resolved = [:]
        geometry.each { field, value -> resolved[field.toString()] = value }
        validateGeometry(key, resolved)
        if (key == 'PIPseq') {
            resolved.cbLen = intParam(params.pipseqCbLength, resolved.cbLen, 'pipseqCbLength')
            resolved.umiLen = intParam(params.pipseqUmiLength, resolved.umiLen, 'pipseqUmiLength')
            resolved.umiStart = resolved.cbStart + resolved.cbLen
        }
        return resolved
    }

    private static Map optionsOf(String key, Map definitions) {
        def chain = []
        def current = key
        while (current != null) {
            def definition = definitions[current]
            if (definition == null) {
                throw new IllegalArgumentException("chemistry '${key}' inherits unknown chemistry '${current}'")
            }
            chain.add(0, current)
            current = definition.inherits?.toString()
            if (chain.size() > definitions.size()) {
                throw new IllegalArgumentException("chemistry '${key}' has a circular inheritance chain")
            }
        }
        def options = [:]
        chain.each { name ->
            def source = definitions[name].starOptions
            if (source != null) {
                source.each { option, value -> options[option.toString()] = value }
            }
        }
        return options
    }

    private static void syncGeometry(Map options, Map geometry) {
        if (geometry == null) {
            return
        }
        GEOMETRY_OPTIONS.each { option, field ->
            if (options.containsKey(option)) {
                options[option] = geometry[field]
            }
        }
    }

    private static String legacyOverride(String key, Map params, Map geometry) {
        def paramName = LEGACY[key]
        if (paramName == null) {
            return null
        }
        def override = params[paramName]
        if (override == null) {
            return null
        }
        if (override instanceof Closure) {
            return override.call(geometry?.cbLen, geometry?.umiLen)?.toString()
        }
        return override.toString()
    }

    private static Map mergeGeometry(String key, Map geometry, Map options) {
        if (geometry == null) {
            return null
        }
        def merged = new LinkedHashMap(geometry)
        GEOMETRY_OPTIONS.each { option, field ->
            def value = options[option]
            if (value != null) {
                merged[field] = value
            }
        }
        validateGeometry(key, merged)
        return merged
    }

    private static void validateGeometry(String key, Map geometry) {
        ['cbStart', 'cbLen', 'umiStart', 'umiLen'].each { field ->
            def value = geometry[field]
            if (value == null) {
                throw new IllegalArgumentException("chemistry '${key}' geometry is missing '${field}'")
            }
            def number
            try {
                number = value.toString().trim().toInteger()
            } catch (NumberFormatException ignored) {
                throw new IllegalArgumentException("chemistry '${key}' geometry '${field}' must be an integer, got '${value}'")
            }
            if (number < 1) {
                throw new IllegalArgumentException("chemistry '${key}' geometry '${field}' must be at least 1, got ${number}")
            }
            geometry[field] = number
        }
        if (geometry.barcodeRead != null) {
            def mate
            try {
                mate = geometry.barcodeRead.toString().trim().toInteger()
            } catch (NumberFormatException ignored) {
                throw new IllegalArgumentException("chemistry '${key}' geometry 'barcodeRead' must be 1 or 2, got '${geometry.barcodeRead}'")
            }
            if (!(mate in [1, 2])) {
                throw new IllegalArgumentException("chemistry '${key}' geometry 'barcodeRead' must be 1 or 2, got ${mate}")
            }
            geometry.barcodeRead = mate
        }
    }

    private static int intParam(Object value, Object fallback, String name) {
        def raw = (value == null) ? fallback : value
        if (raw == null) {
            throw new IllegalArgumentException("--${name} is not set")
        }
        def number
        try {
            number = raw.toString().trim().toInteger()
        } catch (NumberFormatException ignored) {
            throw new IllegalArgumentException("--${name} must be an integer, got '${raw}'")
        }
        if (number < 1) {
            throw new IllegalArgumentException("--${name} must be positive, got ${number}")
        }
        return number
    }
}
