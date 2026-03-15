// Converted from C++ to Rust

use crate::io::variant::vcf_header::VcfHeader;
use crate::io::reference::reference_genome::ReferenceGenome;

pub struct VcfHeaderFactory;

impl VcfHeaderFactory {
    pub fn make_header(
        reference: &ReferenceGenome,
        samples: &[String],
        caller_name: &str,
    ) -> VcfHeader {
        let mut header = VcfHeader::new("VCFv4.3");

        for contig_name in reference.contig_names() {
            let length = reference.contig_size(contig_name);
            header.add_contig(contig_name.clone(), length.map(|s| s as u64));
        }

        for sample in samples {
            header.add_sample(sample.clone());
        }

        header.add_filter("PASS", "All filters passed");
        header.add_filter("LowQual", "Low genotype quality");
        header.add_filter("LowDepth", "Low read depth");

        header.add_info("DP", "1", "Integer", "Total read depth at the locus");
        header.add_info("AF", "A", "Float", "Allele frequency");
        header.add_info("MQ", "1", "Float", "RMS Mapping Quality");
        header.add_info("SOMATIC", "0", "Flag", "Indicates if record is a somatic mutation");

        header.add_format("GT", "1", "String", "Genotype");
        header.add_format("DP", "1", "Integer", "Read depth");
        header.add_format("AD", "R", "Integer", "Allelic depths");
        header.add_format("GQ", "1", "Integer", "Genotype quality");
        header.add_format("PL", "G", "Integer", "Phred-scaled genotype likelihoods");
        header.add_format("PS", "1", "Integer", "Phase set identifier");

        header
    }
}
