use crate::basics::contig_region::ContigRegion;

pub trait Mappable {
    fn mapped_begin(&self) -> u32;
    fn mapped_end(&self) -> u32;
    fn contig_region(&self) -> ContigRegion {
        ContigRegion::new(self.mapped_begin(), self.mapped_end()).unwrap()
    }
}
