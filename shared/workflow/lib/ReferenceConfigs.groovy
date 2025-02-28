import org.yaml.snakeyaml.Yaml

class ReferenceConfigs {
    Map<String, Object> config

    ReferenceConfigs(String yamlFilePath) {
        Yaml yaml = new Yaml()
        this.config = yaml.load(new File(yamlFilePath).text)
    }

    List<Map<String, Object>> getRefs() {
        def validLabels = config.keySet()
        def result = []
        validLabels.each { label ->
            def segments = config[label]?.segments
            segments?.each { segment, details ->
                result << [
                    label  : label,
                    segment: segment,
                    refs   : details.refs,
                    seqs   : details.seqs
                ]
            }
        }
        return result
    }
}
