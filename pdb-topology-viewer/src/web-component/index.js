class PdbTopologyViewer extends HTMLElement {

  static get observedAttributes() {
    return ['entry-id', 'entity-id', 'filter-range', 'chain-id', 'display-style', 'error-style', 'menu-style', 'subscribe-events', 'pvapi'];
  }

  constructor() {
    super();
  }

  validateParams() {
    if(typeof this.entryId == 'undefined' || typeof this.entityId == 'undefined' || this.entryId == null || this.entityId == null) return false;
    return true
  }

  invokePlugin(){
    let paramValidatity = this.validateParams();
    if(!paramValidatity) return

    // create an instance of the plugin
    if(typeof this.pluginInstance == 'undefined') this.pluginInstance = new PdbTopologyViewerPlugin();
    
    let options = {
      entryId: this.entryId,
      entityId: this.entityId,
      entropyId: this.entropyId,
    }

    if(typeof this.chainId !== 'undefined' && this.chainId !== null) options['chainId'] = this.chainId;
    if(typeof this.displayStyle !== 'undefined' && this.displayStyle !== null) options['displayStyle'] = this.displayStyle;
    if(typeof this.errorStyle !== 'undefined' && this.errorStyle !== null) options['errorStyle'] = this.errorStyle;
    if(typeof this.menuStyle !== 'undefined' && this.menuStyle !== null) options['menuStyle'] = this.menuStyle;
    if(typeof this.subscribeEvents !== 'undefined' && this.subscribeEvents !== null) options['subscribeEvents'] = this.subscribeEvents;
    if(typeof this.pvAPI !== 'undefined' && this.pvAPI !== null) options['pvAPI'] = this.pvAPI;
    if(typeof this.filterRange !== 'undefined' && this.filterRange !== null) options['filterRange'] = this.filterRange;
    this.pluginInstance.render(this, options);

  }

  attributeChangedCallback() {
    this.entryId = this.getAttribute("entry-id");
    this.entityId = this.getAttribute("entity-id");
    this.chainId = this.getAttribute("chain-id");
    this.filterRange = this.getAttribute("filter-range");
    this.displayStyle = this.getAttribute("display-style");
    this.errorStyle = this.getAttribute("error-style");
    this.menuStyle = this.getAttribute("menu-style");
    this.subscribeEvents = this.getAttribute("subscribe-events");
    this.pvAPI = (/true/i).test(this.getAttribute("pvapi"));
    this.invokePlugin();
  }

}

export default PdbTopologyViewer;

customElements.define('pdb-topology-viewer', PdbTopologyViewer);