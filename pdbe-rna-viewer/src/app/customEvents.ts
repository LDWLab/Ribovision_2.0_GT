export namespace CustomEvents {

    export function create(eventTypeArr: string[]) {
        /*Create new custom events
        @param eventTypeArr {string[]}: a string array of custom event types
        */
        let eventObj = {} as any;
        for (let ei = 0, el = eventTypeArr.length; ei < el; ei++) {
            let eventType = eventTypeArr[ei];
            let event;
            if (typeof MouseEvent == 'function') {
                // current standard
                event = new MouseEvent(eventType, { 'view': window, 'bubbles': true, 'cancelable': true });
            } else if (typeof document.createEvent == 'function') {
                // older standard
                event = document.createEvent('MouseEvents');
                event.initEvent(eventType, true /* bubbles */, true /* cancelable */);

            }
            eventObj[eventType] = event;
        };
        return eventObj;
    }

    export function dispatchCustomEvent(event: any, eventData: any, targetElement: HTMLElement) {
        /*Dispatch custom event
        @param event {any}: the custom event
        @param eventData {any}: the event data
        @param targetElement {HTMLElement}: the event target
        */
        event['eventData'] = eventData;
        targetElement.dispatchEvent(event);
    }

    export function subscribeToComponentEvents(pluginCtx: any) {
        /* Custom event tracking for other PDB viewers
        @param pluginCtx {any}: the RNA viewer plugin instance
        */

        //Add event listeners for PDB Mol* viewer
        document.addEventListener("PDB.molstar.mouseover", ((e: any) => {
            if (e.eventData && e.eventData.auth_seq_id && e.eventData.auth_asym_id === pluginCtx.options.chainId) {
                pluginCtx.selectResidue(e.eventData.auth_seq_id)
            }
        })),
        document.addEventListener("PDB.molstar.mouseout", ((e: any) => {
            pluginCtx.clearSelection(e.eventData.residueNumber)
        }))
    }

}