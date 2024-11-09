package org.cbioportal.service.impl;
import java.math.BigDecimal;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.apache.commons.lang3.BooleanUtils;
import org.cbioportal.model.MolecularProfile.MolecularAlterationType;
import org.cbioportal.model.meta.GenericAssayMeta;
import org.cbioportal.model.EnrichmentType;
import org.cbioportal.model.GenericAssayBinaryEnrichment;
import org.cbioportal.model.GenericAssayCategoricalEnrichment;
import org.cbioportal.model.GenericAssayMolecularAlteration;
import org.cbioportal.model.MolecularProfileCaseIdentifier;
import org.cbioportal.model.MolecularProfile;
import org.cbioportal.model.Sample;
import org.cbioportal.model.GenericAssayEnrichment;
import org.cbioportal.model.GenomicEnrichment;
import org.cbioportal.model.Gene;
import org.cbioportal.model.GeneMolecularAlteration;
import org.cbioportal.persistence.MolecularDataRepository;
import org.cbioportal.service.ExpressionEnrichmentService;
import org.cbioportal.service.GeneService;
import org.cbioportal.service.GenericAssayService;
import org.cbioportal.service.MolecularProfileService;
import org.cbioportal.service.SampleService;
import org.cbioportal.service.exception.MolecularProfileNotFoundException;
import org.cbioportal.service.util.ExpressionEnrichmentUtil;
import org.cbioportal.service.util.FisherExactTestCalculator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;

@Service
public class ExpressionEnrichmentServiceImpl implements ExpressionEnrichmentService {

    @Autowired
    private MolecularProfileService molecularProfileService;
    @Autowired
    private MolecularDataRepository molecularDataRepository;
    @Autowired
    private GeneService geneService;
    @Autowired
    private ExpressionEnrichmentUtil expressionEnrichmentUtil;
    @Autowired
    private GenericAssayService genericAssayService;
    @Autowired
    private SampleService sampleService;
    @Autowired
    private FisherExactTestCalculator fisherExactTestCalculator = new FisherExactTestCalculator();
    @Override
    // transaction needs to be setup here in order to return Iterable from
    // molecularDataService in fetchCoExpressions
    @Transactional(readOnly = true)
    public List<GenomicEnrichment> getGenomicEnrichments(String molecularProfileId,
                                                         Map<String, List<MolecularProfileCaseIdentifier>> molecularProfileCaseSets, EnrichmentType enrichmentType)
        throws MolecularProfileNotFoundException {
        MolecularProfile molecularProfile = molecularProfileService.getMolecularProfile(molecularProfileId);
        List<MolecularAlterationType> validGenomicMolecularAlterationTypes = Arrays.asList(
            MolecularAlterationType.MICRO_RNA_EXPRESSION, MolecularAlterationType.MRNA_EXPRESSION,
            MolecularAlterationType.MRNA_EXPRESSION_NORMALS, MolecularAlterationType.RNA_EXPRESSION,
            MolecularAlterationType.METHYLATION, MolecularAlterationType.METHYLATION_BINARY,
            MolecularAlterationType.PHOSPHORYLATION, MolecularAlterationType.PROTEIN_LEVEL,
            MolecularAlterationType.PROTEIN_ARRAY_PROTEIN_LEVEL,
            MolecularAlterationType.PROTEIN_ARRAY_PHOSPHORYLATION);
        validateMolecularProfile(molecularProfile, validGenomicMolecularAlterationTypes);
        Iterable<GeneMolecularAlteration> maItr = molecularDataRepository
            .getGeneMolecularAlterationsIterableFast(molecularProfile.getStableId());
        List<GenomicEnrichment> expressionEnrichments = expressionEnrichmentUtil.getEnrichments(molecularProfile,
            molecularProfileCaseSets, enrichmentType, maItr);
        List<Integer> entrezGeneIds = expressionEnrichments.stream().map(GenomicEnrichment::getEntrezGeneId)
            .collect(Collectors.toList());
        Map<Integer, List<Gene>> geneMapByEntrezId = geneService
            .fetchGenes(entrezGeneIds.stream().map(Object::toString).collect(Collectors.toList()), "ENTREZ_GENE_ID",
                "SUMMARY")
            .stream().collect(Collectors.groupingBy(Gene::getEntrezGeneId));
        return expressionEnrichments.stream()
            // filter Enrichments having no gene reference object(this
            // happens when multiple
            // entrez ids map to same hugo gene symbol)
            .filter(expressionEnrichment -> geneMapByEntrezId.containsKey(expressionEnrichment.getEntrezGeneId()))
            .map(expressionEnrichment -> {
                Gene gene = geneMapByEntrezId.get(expressionEnrichment.getEntrezGeneId()).get(0);
                expressionEnrichment.setHugoGeneSymbol(gene.getHugoGeneSymbol());
                return expressionEnrichment;
            }).collect(Collectors.toList());
    }
    @Override
    // transaction needs to be setup here in order to return Iterable from
    // molecularDataRepository in getGenericAssayMolecularAlterationsIterable
    @Transactional(readOnly = true)
    public List<GenericAssayEnrichment> getGenericAssayNumericalEnrichments(String molecularProfileId,
                                                                   Map<String, List<MolecularProfileCaseIdentifier>> molecularProfileCaseSets, EnrichmentType enrichmentType)
        throws MolecularProfileNotFoundException {
        MolecularProfile molecularProfile = molecularProfileService.getMolecularProfile(molecularProfileId);
        validateMolecularProfile(molecularProfile, Arrays.asList(MolecularAlterationType.GENERIC_ASSAY), "SUMMARY");

        Iterable<GenericAssayMolecularAlteration> maItr = molecularDataRepository
            .getGenericAssayMolecularAlterationsIterable(molecularProfile.getStableId(), null, "SUMMARY");
        Map<String, List<MolecularProfileCaseIdentifier>> filteredMolecularProfileCaseSets;
        if (BooleanUtils.isTrue(molecularProfile.getPatientLevel())) {
            // Build sampleIdToPatientIdMap to quick find if a sample has shared patientId with other samples
            List<String> sampleIds = molecularProfileCaseSets.values().stream().flatMap(Collection::stream).map(MolecularProfileCaseIdentifier::getCaseId).collect(Collectors.toList());
            List<String> studyIds = Collections.nCopies(sampleIds.size(), molecularProfile.getCancerStudyIdentifier());
            List<Sample> samples = sampleService.fetchSamples(studyIds, sampleIds, "SUMMARY");
            Map<String, Integer> sampleIdToPatientIdMap = getSampleIdToPatientIdMap(molecularProfile.getCancerStudyIdentifier(), sampleIds);
            // Build filteredMolecularProfileCaseSets
            filteredMolecularProfileCaseSets = new HashMap<>();
            filteredMolecularProfileCaseSets = buildFilteredMolecularProfileCaseSets(molecularProfileCaseSets, sampleIdToPatientIdMap);
        } else {
            filteredMolecularProfileCaseSets = molecularProfileCaseSets;
        }
        List<GenericAssayEnrichment> genericAssayEnrichments = expressionEnrichmentUtil.getEnrichments(molecularProfile,
            filteredMolecularProfileCaseSets, enrichmentType, maItr);
        List<String> getGenericAssayStableIds = genericAssayEnrichments.stream()
            .map(GenericAssayEnrichment::getStableId).collect(Collectors.toList());
        Map<String, GenericAssayMeta> genericAssayMetaByStableId = fetchGenericAssayMeta(getGenericAssayStableIds, molecularProfileId);
        return genericAssayEnrichments.stream().map(enrichmentDatum -> {
            enrichmentDatum.setGenericEntityMetaProperties(
                genericAssayMetaByStableId.get(enrichmentDatum.getStableId()).getGenericEntityMetaProperties());
            return enrichmentDatum;
        }).collect(Collectors.toList());
    }

    @Override
    @Transactional(readOnly = true)
    public List<GenericAssayBinaryEnrichment> getGenericAssayBinaryEnrichments(
        String molecularProfileId,
        Map<String, List<MolecularProfileCaseIdentifier>> molecularProfileCaseSets, EnrichmentType enrichmentType)
        throws MolecularProfileNotFoundException {

        // Validate and fetch molecular profile
        MolecularProfile molecularProfile = getAndValidateMolecularProfile(molecularProfileId, "BINARY");

        // Get the molecular alterations for the provided profile
        Iterable<GenericAssayMolecularAlteration> maItr = molecularDataRepository
            .getGenericAssayMolecularAlterationsIterable(molecularProfile.getStableId(), null, "SUMMARY");

        // Filter the case sets based on molecular profile
        Map<String, List<MolecularProfileCaseIdentifier>> filteredMolecularProfileCaseSets = filterMolecularProfileCaseSets(molecularProfile, molecularProfileCaseSets);

        // Obtain binary enrichments from the utility
        List<GenericAssayBinaryEnrichment> genericAssayBinaryEnrichments = expressionEnrichmentUtil.getGenericAssayBinaryEnrichments(molecularProfile,
            filteredMolecularProfileCaseSets, enrichmentType, maItr);

        // Calculate q-values for enrichments
        calcQValues(genericAssayBinaryEnrichments);

        // Extract stable IDs from binary enrichments
        List<String> getGenericAssayStableIds = genericAssayBinaryEnrichments.stream()
            .map(GenericAssayEnrichment::getStableId).collect(Collectors.toList());

        // Fetch metadata of generic assays by their stable IDs
        Map<String, GenericAssayMeta> genericAssayMetaByStableId = fetchGenericAssayMeta(getGenericAssayStableIds, molecularProfileId);

        // Assign meta properties to each enrichment
        return genericAssayBinaryEnrichments.stream().map(enrichmentDatum -> {
            enrichmentDatum.setGenericEntityMetaProperties(
                genericAssayMetaByStableId.get(enrichmentDatum.getStableId()).getGenericEntityMetaProperties());
            return enrichmentDatum;
        }).collect(Collectors.toList());
    }

    @Override
    @Transactional(readOnly = true)
    public List<GenericAssayCategoricalEnrichment> getGenericAssayCategoricalEnrichments(String molecularProfileId,
                                                                                         Map<String, List<MolecularProfileCaseIdentifier>> molecularProfileCaseSets, EnrichmentType enrichmentType)
        throws MolecularProfileNotFoundException {

        MolecularProfile molecularProfile = getAndValidateMolecularProfile(molecularProfileId, "CATEGORICAL");

        Iterable<GenericAssayMolecularAlteration> maItr = molecularDataRepository
            .getGenericAssayMolecularAlterationsIterable(molecularProfile.getStableId(), null, "SUMMARY");

        Map<String, List<MolecularProfileCaseIdentifier>> filteredMolecularProfileCaseSets = filterMolecularProfileCaseSets(molecularProfile, molecularProfileCaseSets);

        List<GenericAssayCategoricalEnrichment> genericAssayCategoricalEnrichments = expressionEnrichmentUtil.getGenericAssayCategoricalEnrichments(molecularProfile,
            filteredMolecularProfileCaseSets, enrichmentType, maItr);

        calcQValues(genericAssayCategoricalEnrichments);

        List<String> getGenericAssayStableIds = genericAssayCategoricalEnrichments.stream()
            .map(GenericAssayEnrichment::getStableId).collect(Collectors.toList());
        Map<String, GenericAssayMeta> genericAssayMetaByStableId = fetchGenericAssayMeta(getGenericAssayStableIds, molecularProfileId);

        return genericAssayCategoricalEnrichments.stream().map(enrichmentDatum -> {
            enrichmentDatum.setGenericEntityMetaProperties(
                genericAssayMetaByStableId.get(enrichmentDatum.getStableId()).getGenericEntityMetaProperties());
            return enrichmentDatum;
        }).collect(Collectors.toList());
    }

    private MolecularProfile getAndValidateMolecularProfile(String molecularProfileId, String dataType) throws MolecularProfileNotFoundException {
        MolecularProfile molecularProfile = molecularProfileService.getMolecularProfile(molecularProfileId);
        validateMolecularProfile(molecularProfile, Arrays.asList(MolecularAlterationType.GENERIC_ASSAY), "SUMMARY");
        return molecularProfile;
    }

    private void validateMolecularProfile(MolecularProfile molecularProfile, 
                                      List<MolecularAlterationType> validTypes, 
                                      String requiredDataType) throws MolecularProfileNotFoundException {

    if (!validTypes.contains(molecularProfile.getMolecularAlterationType())) {
        throw new MolecularProfileNotFoundException("Invalid MolecularAlterationType for profile: " + molecularProfile.getStableId());
    }
    if (molecularProfile.getMolecularAlterationType().equals(MolecularAlterationType.GENERIC_ASSAY) &&
        !molecularProfile.getDatatype().equals(requiredDataType)) {
        throw new MolecularProfileNotFoundException("Invalid data type for profile: " + molecularProfile.getStableId());
    }
}


    private Map<String, List<MolecularProfileCaseIdentifier>> filterMolecularProfileCaseSets(MolecularProfile molecularProfile, Map<String, List<MolecularProfileCaseIdentifier>> molecularProfileCaseSets) {
        if (BooleanUtils.isTrue(molecularProfile.getPatientLevel())) {
            // If patient level, filter duplicates by patient id
            // For now we only support sample level for samples
            List<String> sampleIds = molecularProfileCaseSets.values().stream().flatMap(Collection::stream).map(MolecularProfileCaseIdentifier::getCaseId).collect(Collectors.toList());
            List<String> studyIds = Collections.nCopies(sampleIds.size(), molecularProfile.getCancerStudyIdentifier());
            List<Sample> samples = sampleService.fetchSamples(studyIds, sampleIds, "ID");
            Map<String, Integer> sampleIdToPatientIdMap = getSampleIdToPatientIdMap(molecularProfile.getCancerStudyIdentifier(), sampleIds);

            Map<String, List<MolecularProfileCaseIdentifier>> filteredMolecularProfileCaseSets = new HashMap<>();
            filteredMolecularProfileCaseSets = buildFilteredMolecularProfileCaseSets(molecularProfileCaseSets, sampleIdToPatientIdMap);
            return filteredMolecularProfileCaseSets;
        } else {
            return molecularProfileCaseSets;
        }
    }
    private List<MolecularProfileCaseIdentifier> filterIdentifiersByUniquePatientId(
        List<MolecularProfileCaseIdentifier> identifiers, Map<String, Integer> sampleIdToPatientIdMap) {

    Set<Integer> patientSet = new HashSet<>();
    List<MolecularProfileCaseIdentifier> uniqueIdentifiers = new ArrayList<>();

    for (MolecularProfileCaseIdentifier caseIdentifier : identifiers) {
        Integer patientId = sampleIdToPatientIdMap.get(caseIdentifier.getCaseId());
        if (!patientSet.contains(patientId)) {
            uniqueIdentifiers.add(caseIdentifier);
            patientSet.add(patientId);
        }
    }
    return uniqueIdentifiers;
}
// Refactored method for setting up filtered molecular profile case sets
private Map<String, List<MolecularProfileCaseIdentifier>> setupFilteredMolecularProfileCaseSets(
        Map<String, List<MolecularProfileCaseIdentifier>> molecularProfileCaseSets, 
        Map<String, Integer> sampleIdToPatientIdMap) {

    Map<String, List<MolecularProfileCaseIdentifier>> filteredSets = new HashMap<>();

    for (Map.Entry<String, List<MolecularProfileCaseIdentifier>> pair : molecularProfileCaseSets.entrySet()) {
        List<MolecularProfileCaseIdentifier> uniqueIdentifiers = filterIdentifiersByUniquePatientId(pair.getValue(), sampleIdToPatientIdMap);
        filteredSets.put(pair.getKey(), uniqueIdentifiers);
    }
    return filteredSets;
}
// Method to get unique identifiers for a single entry based on patient IDs
private List<MolecularProfileCaseIdentifier> getUniqueIdentifiersByPatientId(
        List<MolecularProfileCaseIdentifier> identifiers, Map<String, Integer> sampleIdToPatientIdMap) {

    Set<Integer> patientSet = new HashSet<>();
    List<MolecularProfileCaseIdentifier> uniqueIdentifiers = new ArrayList<>();

    for (MolecularProfileCaseIdentifier caseIdentifier : identifiers) {
        Integer patientId = sampleIdToPatientIdMap.get(caseIdentifier.getCaseId());
        if (patientId != null && !patientSet.contains(patientId)) {
            uniqueIdentifiers.add(caseIdentifier);
            patientSet.add(patientId);
        }
    }
    return uniqueIdentifiers;
}

// Method to build the entire filtered molecular profile case set map
private Map<String, List<MolecularProfileCaseIdentifier>> buildFilteredMolecularProfileCaseSets(
        Map<String, List<MolecularProfileCaseIdentifier>> molecularProfileCaseSets, 
        Map<String, Integer> sampleIdToPatientIdMap) {

    Map<String, List<MolecularProfileCaseIdentifier>> filteredSets = new HashMap<>();
    molecularProfileCaseSets.forEach((key, identifiers) -> 
        filteredSets.put(key, getUniqueIdentifiersByPatientId(identifiers, sampleIdToPatientIdMap))
    );
    return filteredSets;
}
private Map<String, Integer> getSampleIdToPatientIdMap(String studyId, List<String> sampleIds) {
    List<String> studyIds = Collections.nCopies(sampleIds.size(), studyId);
    List<Sample> samples = sampleService.fetchSamples(studyIds, sampleIds, "SUMMARY");
    return samples.stream().collect(Collectors.toMap(Sample::getStableId, Sample::getPatientId));
}


    private Map<String, GenericAssayMeta> getGenericAssayMetaByStableId(List<String> stableIds, String molecularProfileId) {
        return genericAssayService.getGenericAssayMetaByStableIdsAndMolecularIds(stableIds, stableIds.stream().map(sid -> molecularProfileId)
                .collect(Collectors.toList()), "SUMMARY").stream()
            .collect(Collectors.toMap(GenericAssayMeta::getStableId, Function.identity()));
    }
    private Map<String, GenericAssayMeta> fetchGenericAssayMeta(List<String> stableIds, String molecularProfileId) {
        return genericAssayService
            .getGenericAssayMetaByStableIdsAndMolecularIds(stableIds, 
                Collections.nCopies(stableIds.size(), molecularProfileId), "SUMMARY")
            .stream()
            .collect(Collectors.toMap(GenericAssayMeta::getStableId, Function.identity()));
    }
    
    private <T extends GenericAssayEnrichment> void calcQValues(List<T> enrichments) {
        // Sort enrichments by pValue
        Collections.sort(enrichments, GenericAssayEnrichment::compare);
        BigDecimal[] pValues = enrichments.stream().map(T::getpValue).toArray(BigDecimal[]::new);
        BigDecimal[] qValues = fisherExactTestCalculator.calcqValue(pValues);
        // Assign q-values to enrichments
        for (int i = 0; i < enrichments.size(); i++) {
            enrichments.get(i).setqValue(qValues[i]);
        }
    }
    
    private void validateMolecularProfile(MolecularProfile molecularProfile,
                                          List<MolecularAlterationType> validMolecularAlterationTypes) throws MolecularProfileNotFoundException {
        if (!validMolecularAlterationTypes.contains(molecularProfile.getMolecularAlterationType())) {
            throw new MolecularProfileNotFoundException(molecularProfile.getStableId());
        }
    }
}
